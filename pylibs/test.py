from pyIPSIG import *
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

import matplotlib.pyplot as plt

import time

def getRandomE(N, Cr=0.15, show = False):
	e = None
	goon = True
	
	while goon:
		e = np.zeros((N,N), dtype = bool)
		for i in range(N):
			for j in range(i+1,N):
				if i == 0 and j == N-1:		## otherwise it would break the minpath algorithm
					continue
				r = np.random.rand()
				if r < Cr:
					e[i,j] = 1
			
		elist = []
		for i in range(N):
			for j in range(N):
				if e[i,j]:
					elist.append([i,j])
						
		g = nx.Graph()
		g.add_nodes_from([x for x in range(N)])
		g.add_edges_from(elist)
		
		if not nx.is_connected(g):
			continue
		
		paths = dict(nx.shortest_path_length(g))
		maxind = [None, None]; maxd = 0
		for x in paths:
			for y in paths[x]:
				d = paths[x][y]
				if d > maxd:
					maxd = d
					maxind = [x,y]
#		print('maxind')
#		print(maxind)
		if maxd>0:
			a = maxind[0]; b = maxind[1]
			hold = e[0,:].copy()
			e[0,:] = e[a,:].copy()
			e[a,:] = hold.copy()
			
			hold = e[:,b].copy()
			e[:,b] = e[:,N-1].copy()
			e[:,N-1] = hold.copy()
			
		elist = []
		for i in range(N):
			for j in range(N):
				if e[i,j]:
					elist.append([i,j])
			
						
		g = nx.Graph()
		g.add_nodes_from([x for x in range(N)])
		g.add_edges_from(elist)
		
		if not nx.is_connected(g):
			continue



		out = []
		for i in range(N):
			out.append([])
		for e in elist:
			out[e[0]].append(e[1])
			out[e[1]].append(e[0])

		if show:			
			nx.draw(g, labels = {x:str(x) for x in range(N)})
			plt.show(0)
			
			return out, g
		
		
		break
		
	return out

N = 25

times = []
Crs = [0.1, 0.2, 0.3]
Ns = [15,25,30,35,40]
for Cr in Crs:
	for N in Ns:
		for i in range(50):
			print([N,Cr,i])
			s1 = time.time()
			e = getRandomE(N+2, Cr)
			g = Graph(e)
#			print('building graph, time %.3f' % (time.time()-s1))
			
			
			start = time.time()	
			s1 = time.time()
			s = System(g, [0]*N)
			#s = System(g, np.random.randint(4, size = N))
#			print('init, time %.3f' % (time.time()-s1))
			
			i = 1
			while s.branch:
				s1 = time.time()
				s.iterate()
#				print('\titerace %d, b: %d, time %.3f' % (i,len(s.branch), time.time()-s1))
				i+=1
				#print(s.sig0)
				#print(s.sig1)
				#print('------------------------')
			end = time.time()
			times.append([end - start, N, Cr])
	
	
import seaborn as sns
import pandas as pd
data = pd.DataFrame(times, columns = ['time','system size', 'Pedge'])
data['log time'] = np.log10(data['time'])
sns.violinplot(x='system size',y='log time', hue='Pedge',data=data, inner = 'points')
plt.show(0)
