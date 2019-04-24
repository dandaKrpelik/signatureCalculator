#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include "omp.h"
#include "pChrono.h"
#include "graph.h"

#include <cstdint>

using namespace std;

using posvec = vector<size_t>;
//~ struct bigInt	
//~ {
	//~ uint64_t bits[4];	//256bin unisgned int
//~ };
using bigInt = __uint128_t;
using indexType = size_t;
using dick = map<indexType, bigInt>;
using fdick = map<indexType, float>;


size_t binomialCoeff(size_t n, size_t k) 
{ 
  // Base Cases 
  if (k==0 || k==n) 
    return 1; 
  
  // Recur 
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k); 
} 

void addToDick(dick* d, indexType pos, bigInt val)
{
	dick::iterator it = d->find(pos);
	if (it == d->end())
		d->insert(make_pair(pos, val));
	else
	{
		bigInt newval = (*d)[pos] + val;
		(*d)[pos] = newval;
	}	
}



void printDick(dick d)
{
	dick::iterator it = d.begin();
	while(it != d.end())
	{
		cout<< "("<<it->first << " : " << it->second << "), ";
		it++;
	}
	cout << endl;
}


class PositionIterator
// to iterate through tensor
{
	private:
	posvec lim;
	posvec val;
	size_t n;
	
	public:
	bool finished;
	PositionIterator(posvec limin)
	{
		lim = limin;
		n = lim.size();
		val.resize(n, 0);
		finished = false;
	}
	void operator++()
	{
		finished = true;
		for (size_t i = 0; i < n; i++)
		{
			val[i] += 1;
			if (val[i] <= lim[i])
			{
				finished = false;
				break;
			}
			else
			{
				val[i] = 0;
			}
		}
	}
	posvec operator*()
	{
		return val;
	}
	void copypos(posvec* pos)
	{
		copy(val.begin(), val.end(), pos->begin());
	}
};


class Branch
{
	public:
	vector<bool> mask;		//which components are not failed
	posvec ones;			//which components need to remain functional due to the state space decomposition
	posvec minpath;
	
	Branch(vector<bool> maskin, posvec onesin, posvec minpathin)
	{
		mask = maskin;
		ones = onesin;
		minpath = minpathin;
	}
};


class System
{
	public:
	
	Graph* g;	//RBD
	size_t N;	//no. components
	size_t K;	//no. types
	posvec M;	//amount of components of resp. types
	posvec C;	//component types
	
	posvec mul;	//multiplication coefficient for calculating position in multi-dimensional array
	
	vector<Branch*> branch;
	
	vector<indexType> sig1;		//counting ones and zeros
	vector<indexType> sig0;
		
	vector<posvec> binCo;		//precalculate binomial coefficients
	
	
	bool validSig = 0;
	fdick LowSig;
	fdick HiSig;
		
	~System()
	{
//		delete g;
		for (Branch* b : branch)
			delete b;
	}
		
	System(Graph* gin, posvec Cin)
	{
		g = gin;
		C = Cin;
		
		// from type assigments C, create vector of number of components of various types, M
		N = C.size();
		
		dick typeCounts;
		size_t maxType = 0;
		for (size_t i = 0 ; i < N ; i++)
		{
			if (C[i] > maxType)	
				maxType = C[i];
			addToDick(&typeCounts, C[i], 1);
		}
		//cout << maxType << endl;
		//printDick(typeCounts);
				
		K = maxType+1;	// number of types
		M.resize(K,0);
		for (size_t i = 0; i< maxType+1; i++)
		{
			dick::iterator it = typeCounts.find(i);
			if (it != typeCounts.end())
				M[i] = typeCounts[i];	
		}
		
		mul.resize(K, 1);
		size_t maxM = 0;
		for (int i = 0 ; i < K; i++)
		{
			if (M[i] > maxM) maxM = M[i];
			if (i==0) continue;
			mul[i] = (M[i-1]+1)*mul[i-1];
		}
		
		//precalculate bin.coefs
		maxM += 1;
		for (size_t i = 0 ; i < maxM; i++)
		{
			posvec fnc(i+1, 1);
			if (i > 0)
				for (size_t j = 1; j < i; j++)
				{
					fnc[j] = binCo[i-1][j-1] + binCo[i-1][j];	// so its a 2d array s.t. BC[i,j] = binCoef(i,j)
				}
			binCo.push_back(fnc);
		}
		
		//signature size
		size_t sig_n = mul[K-1]*(M[K-1]+1);
		sig1.resize(sig_n,0);
		sig0.resize(sig_n,0);
		
		branch.resize(1);
		posvec minpath = minPath_plain(*g);
		vector<bool> firstMask(N, true);
		posvec firstOnes;
		branch[0] = new Branch(firstMask, firstOnes, minpath);
		
		registerBranches();
	}
	
	bool eval(vector<bool> x)
	{
		posvec minpath = minPath(*g, x);
		return !minpath.empty();
	}
	
	vector<float> SIGp(posvec pos)
	{
		indexType index = pos2index(pos);
		return SIGi(index);
	}
	vector<float> SIGi(indexType index)
	{
		if (!validSig) calcSIG();
		
		vector<float> P(2);
		fdick::iterator it = LowSig.find(index);
		if (it == LowSig.end())
		{
			P[0] = 0.;
		}
		else
		{
			P[0] = it->second;
		}
		
		it = HiSig.find(index);
		if (it == LowSig.end())
		{
			P[1] = 1.;
		}
		else
		{
			P[1] = it->second;
		}
		return P;
	}
	
	void calcSIG()
	{
		if (validSig) return;
		
		LowSig.clear();
		HiSig.clear();
		
		PositionIterator it(M);
		posvec P(K);
		indexType index;
		while (!it.finished)
		{
			it.copypos(&P);
			index = pos2index(P);
			bigInt Nsum = 1;
			float Plow = sig1[index];
			float Phi = sig0[index];
			for (size_t ki = 0; ki < K; ki++)
			{
				//Nsum *= binCo[M[ki]][pos[ki]];
				Plow /= binCo[M[ki]][P[ki]];
				Phi /= binCo[M[ki]][P[ki]];
			}			
			Phi = 1. - Phi;
			if (Plow > 0) LowSig.insert(make_pair(index, Plow));
			if (Phi > 0 ) HiSig.insert(make_pair(index, Phi));
			
			++it;
		}
		validSig = 1;
	}

	indexType pos2index(posvec pos)
	{
		indexType out = 0;
		for (int i = 0; i < K; i++)
		{
			out += mul[i] * pos[i];
		}
		
		return out;
	}
	
	posvec index2pos(indexType index)
	{
		posvec out(K,0);
		bigInt iter = index;
		for (int i = K-1; i > -1; i--)
		{
			out[i] = iter/mul[i];
			iter -= out[i]*mul[i];
		}
		
		return out;
	}
	
	void addFreeCombinations(dick* target, Branch* b)
	// when we arrive to certain failed of functional states -> how many states are there?
	{
		vector<bool> knownOnes(N, false);
		posvec freeM = M;
		size_t pos = 0;
		
		size_t ones_n = b->ones.size();
		for (size_t i = 0; i < ones_n ; i++)
		{
			size_t label = b->ones[i];
			knownOnes[ label ] = true;
			pos += mul[ C[label] ];
		}
		
		size_t path_n = b->minpath.size();
		for (size_t i = 0 ; i < path_n; i++)
		{
			size_t label = b->minpath[i];
			if (!knownOnes[label])
				pos += mul[ C[label] ];
			knownOnes[ label ] = true;
		}

		for (size_t i = 0; i< N; i++)
		{
			if (knownOnes[i] || !b->mask[i])
				freeM[ C[i] ] -= 1;
		}
		
		PositionIterator it(freeM);
		posvec P(K);
		while (!it.finished)
		{
			it.copypos(&P);
			bigInt dpos = 0;
			bigInt val = 1;
			for (size_t i = 0; i < K; i++)
			{
				dpos += mul[i]*P[i];
				//val *= binomialCoeff(freeM[i], P[i]);
				val *= binCo[ freeM[i] ][ P[i] ];
			}
			addToDick(target, pos+dpos, val);
			++it;
		}
	}
	
	int registerBranches()
	// cycles trought all the branches in `branch' to add certain states to sig1
	{
		size_t n = branch.size();
		#pragma omp parallel
		{
			dick t_sig1;		//thread ones counter
			
			#pragma omp for schedule(dynamic, 1)
			for (size_t iit = 0; iit < n ; iit++)
			{
				Branch* b = branch[iit];								
				//and add signature contribution for minpath
				addFreeCombinations(&t_sig1, b);
				
			}
			#pragma omp critical
			{
				dick::iterator it = t_sig1.begin();
				while(it != t_sig1.end())
				{
					//cout << "adding, " << it->first <<" with " << it->second<<endl;
					sig1[it->first] += it->second;
					it++;
				}
				
			}
		}
		return 0;
	}
	
	int iterate()
	{
		if (branch.size() == 0)
			return 0;
		
		validSig = 0;
		if (branch.size() == 1)
			return iterate1();
		else
			return iterateM();
	}
	
	int iterateM()
	{
		size_t n = branch.size();
		vector<Branch*> allNewBranches;
		#pragma omp parallel
		{
			vector<Branch*> newbranch;
			//dick t_sig1;		//thread ones and zero counters
			dick t_sig0;
			
			#pragma omp for schedule(dynamic)
			for (size_t iit = 0; iit < n ; iit++)
			{
				Branch* b = branch[iit];
				size_t path_n = b->minpath.size();
				
				vector<bool> mask = b->mask;
				posvec cuts;
				
				vector<posvec> minpaths(path_n);
				posvec ones = b->ones;
								
				for (size_t i = 0; i < path_n; i++)
				{
					size_t label = b->minpath[i];
					
					if (has(ones, label))		//this component creates cutset (has been discovered in previous branches) -> no need to look for minpath without it
						continue;
					
					mask[label] = false;
					posvec subminpath = minPath(*g, mask);
					minpaths[i] = subminpath;
					
					if (subminpath.empty() && !has(ones, label))	//this component creates cutset
						cuts.push_back(label);
						
					mask[label] = true;
					
					
					//cleanup
					subminpath.clear();
				}
				
				//~ //add cutsets to sig0
				size_t cut_n = cuts.size();
				posvec foo(0);	// represents empty path
				for (size_t i = 0; i < cut_n; i++)
				{
					size_t label = cuts[i];
					mask[label] = false;					
					Branch* nB = new Branch(mask, ones, foo);
					addFreeCombinations(&t_sig0, nB);
					delete nB;

					mask[label] = true;
					ones.push_back(label);
				}
				
				ones = b->ones;
				ones.insert(ones.end(), cuts.begin(), cuts.end());		//ensure that the cut sets wont be investigated later to save resources

				for (size_t i = 0; i < path_n; i++)
				//here, we split the state space into 0xxx,10xx,110x,1110 ... 111111xxxxx get treated separatedly in 'registerBranches'
				{
					size_t label = b->minpath[i];
					mask[label] = false;

					if (!has(ones, label))	//if its worth it
						newbranch.push_back(new Branch(mask, ones, minpaths[i]));
					
					mask[label] = true;
					if (!has(ones, label))		//if its not there yet
						ones.push_back(label);
				}
				
			}
			#pragma omp critical
			{
				allNewBranches.insert(allNewBranches.end(), newbranch.begin(), newbranch.end());
				dick::iterator it = t_sig0.begin();
				while(it != t_sig0.end())
				{
					sig0[it->first] += it->second;
					it++;
				}
			}
		}
		
		for (Branch* b : branch)
			delete b;
		branch = allNewBranches;
		registerBranches();
		
		return 0;
	}
	
	int iterate1()
	{
		vector<Branch*> allNewBranches;
		
		Branch* b = branch[0];	//only one branch here
		size_t path_n = b->minpath.size();
		posvec cuts(path_n, -1);
		vector<posvec> minpaths(path_n);
		posvec ones = b->ones;
		
		posvec foo(0);	// represents empty path

		posvec validCuts;
		size_t validCuts_n;
		vector<posvec> ones_all;

		#pragma omp parallel
		{
			//dick t_sig1;		//thread ones and zero counters
			dick t_sig0;
			
			vector<bool> t_mask = b->mask;
			
			#pragma omp for schedule(dynamic)	
			for (size_t i = 0; i < path_n; i++)
			{
				size_t label = b->minpath[i];
					
				if (has(ones, label))		//this component creates cutset (has been discovered in previous branches) -> no need to look for minpath without it
					continue;
					
				t_mask[label] = false;
				posvec subminpath = minPath(*g, t_mask);
				minpaths[i] = subminpath;
					
				if (subminpath.empty() && !has(ones, label))	//this component creates cutset
					cuts[i] = label;
						
				t_mask[label] = true;
			}

			#pragma omp barrier
			#pragma omp single
			{
				for (size_t i = 0; i < path_n; i++)
				{
					if (cuts[i] < N)
					{
						posvec o = b->ones;
						for (int j = 0; j < validCuts.size(); j++)
						{
							o.push_back( validCuts[j] );
						}
						validCuts.push_back(cuts[i]);
						ones_all.push_back(o);
					}
				}
				validCuts_n = validCuts.size();
			}

			#pragma omp barrier

			#pragma omp for schedule(dynamic)
			for (size_t i = 0; i < validCuts_n; i++)
			{
				size_t label = validCuts[i];
				
				t_mask[label] = false;					
				Branch* nB = new Branch(t_mask, ones_all[i], foo);
				addFreeCombinations(&t_sig0, nB);
				delete nB;
				
				t_mask[label] = true;
			}
				
			#pragma omp single
			{
				ones = b->ones;
				ones.insert(ones.end(), validCuts.begin(), validCuts.end());		//ensure that the cut sets wont be investigated later to save resources
				for (size_t i = 0; i < path_n; i++)
				//here, we split the state space into 0xxx,10xx,110x,1110 ... 111111xxxxx get treated separatedly in 'registerBranches'
				{
					size_t label = b->minpath[i];
					t_mask[label] = false;

					if (!has(ones, label))	//if its worth it
						allNewBranches.push_back(new Branch(t_mask, ones, minpaths[i]));
						
					t_mask[label] = true;
					if (!has(ones, label))		//if its not there yet
						ones.push_back(label);
				}
				
			}
			
			#pragma omp critical
			{
				dick::iterator it = t_sig0.begin();
				while(it != t_sig0.end())
				{
					sig0[it->first] += it->second;
					it++;
				}
			}
		}
		
		delete b;
		branch = allNewBranches;
		registerBranches();
	
		return 0;
	}
};
