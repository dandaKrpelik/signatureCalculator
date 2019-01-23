#include "omp.h"
#include "graph.h" 
#include "system.h"
#include "pChrono.h"

#include <iostream>
#include <vector>


using namespace std;



int main()
{
	helloWorld(17);
  
	vector<bool> x = {1,1,0,1};
	vector<bool> y = {1,1,1,1};

	
	if (x[0])
	{
		cout << "x0 is"<<endl;
	}
  
	if (!has(x,false)) cout << "x is pos"  <<endl;
	if (!has(y,false)) cout << "y is pos"  <<endl;
  
	pChrono.start("opm test");
  
	size_t N = 1000000;
	vector<float> out(N);
	#pragma omp parallel 
	{
	cout << omp_get_thread_num() << endl;
	#pragma omp for
		for (size_t i = 0; i<1000*N; i++) out[i%N] = i*1.1;
	}
  
	pChrono.end("omp test end");
  
	posvec M(3);
	M[0] = 3;
	M[1] = 2;
	M[2] = 3;
	
	PositionIterator it(M);
	while (!it.finished)
	{
		printVector(*it);
		++it;
	}
	
	return 0;
}
