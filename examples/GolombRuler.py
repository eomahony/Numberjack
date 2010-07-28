#! /usr/bin/env python

from Numberjack import *

param = input({'solver':'Mistral', 'marks':6})

ruler = (0,1,3,6,11,17,25,34,44,55,72,85,106,127)
m = param['marks']
n = 2**(m-1)

marks = VarArray(m,n,'m')
distance = [marks[i] - marks[j] for i in range(1,m) for j in range(i)]

model = Model(
    Minimise( marks[-1] ), #objective function
    
    [marks[i-1] < marks[i] for i in range(1,m)],
    AllDiff(distance),
    marks[0] == 0, # symmetry breaking

    [distance[i*(i-1)/2+j] >= ruler[i-j] for i in range(1,m) for j in range(0,i-1) if (i-j < m-1)]
)

s = model.load( param['solver'], marks )
s.setHeuristic('Impact')
s.setTimeLimit(10)
s.setVerbosity(2)
s.solve()

print "Search nodes:", s.getNodes(), "Search time:", s.getTime()
print  marks




