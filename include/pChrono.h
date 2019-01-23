#pragma once

#include <iostream>
#include <vector>
#include <map>

using namespace std;

class pclock
{
	vector<clock_t> time;
	vector<ulong> inds;
	
  map<string,pair<float, size_t> > stat;
  
  void reg(string k, float val)
  {
    pair<float, size_t> a(val,1);
    std::pair<std::map<string ,  pair<float, size_t> >::iterator,bool> ret;
    
//    #pragma omp critical
//    {
      pair<string, pair<float, size_t> > n(k,a);
      ret = stat.insert( n );
      if (ret.second==false)
      {
        stat[k].first += val;
        stat[k].second += 1;
      }
//    }
  }
  
  public:
  
  float getStat(string key)
  {
    std::map<string, pair<float, size_t> >::iterator ret;
    ret = stat.find(key);
    
    if(ret  == stat.end() ) return -1;
    return ret->second.first / ret->second.second;
  }
  
	void start(string msg = "", bool silent = false)
  {
    inds.push_back(time.size()); time.push_back(clock());  
    if(!silent)
      cout<<msg<<" "<<endl;
  }
	void tick(string msg = "", bool silent = false)
  {
    time.push_back(clock()); 
    float val = float( *(time.end()-1) - *(time.end()-2) )/CLOCKS_PER_SEC;
    reg(msg, val);
    if(!silent)
      cout << msg <<" "<< val << " s"<<endl;
  }
	void end(string msg = "", bool silent = false)
  {
    time.push_back(clock());
    float val = float( *(time.end()-1) - time[inds.back()] )/CLOCKS_PER_SEC;
    if(!silent)
      cout << msg <<" "<< val << " s"<<endl; 
    inds.pop_back();
  }
	void reset()
  {time.clear();}
}pChrono;


static void chStart(string s)
{
    pChrono.start(s);
}
static void chEnd(string s)
{
    pChrono.end(s);
}
