#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include "omp.h"
#include "pChrono.h"

using namespace std;
using duo = pair<size_t,size_t>;






// auxiliary functions

ostream&
operator<<( std::ostream& dest, __uint128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}




template <typename T>
bool has(vector<T> x, T val)
{
	return ( find(x.begin(), x.end(), val) != x.end() );
}

void printVector(vector<int> v)
{
	int n = v.size();
	for (int i = 0; i < n ; i++)
		cout << v[i] << ", ";
	cout << endl;
}
void printVector(vector<__uint128_t> v)
{
	int n = v.size();
	for (int i = 0; i < n ; i++)
		cout << v[i] << ", ";
	cout << endl;
}
void printVector(vector<size_t> v)
{
	int n = v.size();
	for (int i = 0; i < n ; i++)
		cout << v[i] << ", ";
	cout << endl;
}

void printVector(vector<duo> v)
{
	size_t n = v.size();
	for (size_t i = 0; i < n ; i++)
		cout << "(" << v[i].first <<","<<v[i].second<< ") , ";
	cout << endl;
}

void printVector(vector<bool> v)
{
	int n = v.size();
	for (int i = 0; i < n ; i++)
		cout << ( v[i] ? "T" : "F" ) << ", ";
	cout << endl;
}

class Graph
//graph data structure
{
	public:
	vector<vector<size_t> > E;	//edges
	int N;	//no. of graph vertices (incl. start and terminal); labeling 0,1,....,N-1
	int _N; //no.  --//--   (excl. start and terminal)
	
	
	Graph (vector<vector<size_t> > Ein)
	// 
	{
		N = Ein.size();
		_N = N - 2;	
		
		E.resize(N);
		for (size_t i = 0 ; i < N ; i++)
		{
			E[i] = Ein[i];
		}
	}
};


static vector<size_t> minPath(Graph g, vector<bool> mask)
// returns set of vertices on minimal path; indexes -1 to use in sig calc
// if start and terminal are connected -> seg. error
{	
	size_t N = g.N;
	vector<int> label(N, -1);	// it's flooding algorithm, this vector holds which liquid has gotten into
	vector<size_t> history(N, -1);	// holds from which vertex the current vertex was flooded (for path reconstruction)
	
	label[0] = 0;
	label[N-1] = 1;
	
	vector<duo> start_neighbours, term_neighbours;
	size_t encounter = -1;		//once I'll encounter flooded vertex, whose predecessor is saved in 'history', I save its id here
	size_t otherNeighbour = -1;		//once I'll encounter flooded vertex, whose predecessor is saved in 'history', I save the other one here (for path reconstruction)
	
	
	//init
	size_t n = g.E[0].size();
	for (size_t i = 0; i < n; i++)
	{
		size_t lab = g.E[0][i];
		if (!mask[ lab -1 ]) continue;
		duo tba(lab, 0);
		start_neighbours.push_back( tba );
	}
	n = g.E[N-1].size();
	for (size_t i = 0; i < n; i++)
	{
		size_t lab = g.E[N-1][i];   
		if (!mask[ lab -1 ]) continue;
		duo tba(lab, N-1);
		term_neighbours.push_back( tba );
	}
	
	
	bool done = false;
	vector<duo> newvector;
	//FLOOD
	while(!done)
	{
//DBG		printVector(start_neighbours);
//DBG		printVector(term_neighbours);
		if (start_neighbours.size() == 0 && term_neighbours.size() == 0) break;
		
		newvector.clear();
		n = start_neighbours.size();

		for (size_t i = 0; i < n; i++)
		{
			size_t lab = start_neighbours[i].first;
			size_t prev = start_neighbours[i].second;
			if (label[lab] == 0) continue;		//if its already flooded, skip it
			
			if (label[lab] == 1)
			//break condition - terminal encountered
			{
				encounter = lab;
				otherNeighbour = prev;
				done = true;
				break;
			}
			
			label[lab] = 0;
			history[lab] = prev;
			
			size_t m = g.E[lab].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t nlab = g.E[lab][j];
				if (nlab > 0 && nlab < N-1)
					if (!mask[nlab -1 ]) continue;
				duo tba(nlab, lab);
				newvector.push_back( tba );
			}
		}
		start_neighbours = newvector;
		newvector.clear();
		n = term_neighbours.size();
		for (size_t i = 0; i < n; i++)
		{
			size_t lab = term_neighbours[i].first;   
			size_t prev = term_neighbours[i].second;			
			if (label[lab] == 1) continue;		//if its already flooded, skip it
			
			if (label[lab] == 0)
			//break condition - start encountered
			{
				encounter = lab;
				otherNeighbour = prev;
				done = true;
				break;
			}
			
			label[lab] = 1;
			history[lab] = prev;
			
			size_t m = g.E[lab].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t nlab = g.E[lab][j];
				if (nlab > 0 && nlab < N-1)
					if (!mask[nlab -1 ]) continue;
				duo tba(nlab, lab);
				newvector.push_back( tba );
			}
		}
		term_neighbours = newvector;
	}
	
	
	//RECONSTRUCT
	vector<size_t> path;
	if (!done)		return path;		// there is no path
	
	size_t towardsStart = label[ encounter ] == 0 ? encounter : otherNeighbour;
	size_t towardsTerm = label[ encounter ] == 1 ? encounter : otherNeighbour;

	while(towardsStart > 0)
	{
		path.push_back(towardsStart-1);
		towardsStart = history[towardsStart];
	}
	while(towardsTerm < N-1)
	{
		path.push_back(towardsTerm-1);
		towardsTerm = history[towardsTerm];
	}
		
		
	return path;
	
}

static vector<size_t> minPath_plain(Graph g)
// returns set of vertices on minimal path; indexes -1 to use in sig calc
// if start and terminal are connected -> seg. error
{	
	size_t N = g.N;
	vector<int> label(N, -1);	// it's flooding algorithm, this vector holds which liquid has gotten into
	vector<size_t> history(N, -1);	// holds from which vertex the current vertex was flooded (for path reconstruction)
	
	label[0] = 0;
	label[N-1] = 1;
	
	vector<duo> start_neighbours, term_neighbours;
	size_t encounter = -1;		//once I'll encounter flooded vertex, whose predecessor is saved in 'history', I save its id here
	size_t otherNeighbour = -1;		//once I'll encounter flooded vertex, whose predecessor is saved in 'history', I save the other one here (for path reconstruction)
	
	
	//init
	size_t n = g.E[0].size();
	for (size_t i = 0; i < n; i++)
	{
		size_t lab = g.E[0][i];
		duo tba(lab, 0);
		start_neighbours.push_back( tba );
	}
	n = g.E[N-1].size();
	for (size_t i = 0; i < n; i++)
	{
		size_t lab = g.E[N-1][i];   
		duo tba(lab, N-1);
		term_neighbours.push_back( tba );
	}
	
	
	bool done = false;
	vector<duo> newvector;
	//FLOOD
	while(!done)
	{
//DBG		printVector(start_neighbours);
//DBG		printVector(term_neighbours);
		if (start_neighbours.size() == 0 && term_neighbours.size() == 0) break;
		
		newvector.clear();
		n = start_neighbours.size();

		for (size_t i = 0; i < n; i++)
		{
			size_t lab = start_neighbours[i].first;
			size_t prev = start_neighbours[i].second;
			if (label[lab] == 0) continue;		//if its already flooded, skip it
			
			if (label[lab] == 1)
			//break condition - terminal encountered
			{
				encounter = lab;
				otherNeighbour = prev;
				done = true;
				break;
			}
			
			label[lab] = 0;
			history[lab] = prev;
			
			size_t m = g.E[lab].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t nlab = g.E[lab][j];
				duo tba(nlab, lab);
				newvector.push_back( tba );
			}
		}
		start_neighbours = newvector;
		newvector.clear();
		n = term_neighbours.size();
		for (size_t i = 0; i < n; i++)
		{
			size_t lab = term_neighbours[i].first;   
			size_t prev = term_neighbours[i].second;			
			if (label[lab] == 1) continue;		//if its already flooded, skip it
			
			if (label[lab] == 0)
			//break condition - start encountered
			{
				encounter = lab;
				otherNeighbour = prev;
				done = true;
				break;
			}
			
			label[lab] = 1;
			history[lab] = prev;
			
			size_t m = g.E[lab].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t nlab = g.E[lab][j];
				duo tba(nlab, lab);
				newvector.push_back( tba );
			}
		}
		term_neighbours = newvector;
	}
	
	
	//RECONSTRUCT
	vector<size_t> path;
	if (!done)		return path;		// there is no path
	
	size_t towardsStart = label[ encounter ] == 0 ? encounter : otherNeighbour;
	size_t towardsTerm = label[ encounter ] == 1 ? encounter : otherNeighbour;

	while(towardsStart > 0)
	{
		path.push_back(towardsStart-1);
		towardsStart = history[towardsStart];
	}
	while(towardsTerm < N-1)
	{
		path.push_back(towardsTerm-1);
		towardsTerm = history[towardsTerm];
	}
		
		
	return path;
	
}

static int helloWorld(float a = 1)
{
	cout << a << endl;
	pChrono.start("testing chrono");
	pChrono.end("ending chrono");
	return 0;
};
