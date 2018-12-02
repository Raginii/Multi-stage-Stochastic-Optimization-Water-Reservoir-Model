from gurobipy import *
import pandas as pd

import project_data as data

t = data.t
w = data.w
ft = data.ft 
h = data.h
demand = data.demand
exchange = data.exchange
wlist = data.wlist
wcapacity = data.wcapacity
flist = data.flist
wprices = data.wprices
winflow = data.winflow
swinflow = data.swinflow
resmax = data.resmax
resmin = data.resmin

wt = data.wt
# last = []
# first = []
# p=0
# for i in w:
#     for j in h:
#         wt.append((i,t[p]))
#         if h.index(j)==0:
#             last.append((i,t[p]))
#         if h.index(j)==len(h)-1:
#             first.append((i,t[p]))
#         p = p+1    

# wt = tuplelist(wt)
# first = tuplelist(first)
# last = tuplelist(last)

inflow = data.inflow
capacity = data.capacity
gencost = data.gencost
capacityKeys = capacity.keys()

##########################################################
#Creating Model
m = Model('DeterministicDam')

#Creating variables
gap = m.addVars(t,name = 'gap') ##Gap generation in T
res = m.addVars(t,name = 'res') ##reservoir energy level at end in T
spill = m.addVars(t,name = 'spill') ##Spillafe in hour t
x = m.addVars(t,ft,name = 'x') ##Generated power in hour t with ft
slackup = m.addVars(t,name = 'slackup') ##Slack for the upper level
slacklo = m.addVars(t,name = 'slacklo') ##Slack for the lower level
                   
#update model
m.update()

#Creating constraints
m.addConstrs(((res[i] - res[t[t.index(i)-1]] + x[i,'Hydro'] + spill[i] == inflow[i]) for i in t) , "Hydraulic continuity" )
m.addConstrs(((quicksum(x[i,f] for f in ft if (i,f) in capacityKeys)  + gap[i] >= demand[i] - exchange[i]) for i in t) , "Demand Satisfaction" )
m.addConstrs(((-res[i] + slackup[i] >= -resmax[i]) for i in t) , "Maximum reservoir level" )
m.addConstrs(((res[i] + slacklo[i] >= resmin[i]) for i in t) , "Minimum reservoir level" )
m.addConstrs(((-x[i,f] >= -capacity[i,f]) for i in t for f in ft if (i,f) in capacityKeys) , "Minimum reservoir level" )

#Setting objective
m.setObjective(
        (
        quicksum(gencost[i,f]*x[i,f] for i in t for f in ft if (i,f) in capacityKeys)
        + 1000*quicksum(gap[i] for i in t)
        + 10e6*quicksum(slackup[i]+slacklo[i] for i in t)
        )
        ,GRB.MINIMIZE)

m.update()


#optimize
m.optimize()

status = m.status
if status == GRB.Status.UNBOUNDED:
    print('The model cannot be solved because it is unbounded')
if status == GRB.Status.OPTIMAL:
    print('The optimal objective is %g' % m.objVal)
if status != GRB.Status.INF_OR_UNBD and status != GRB.Status.INFEASIBLE:
    print('Optimization was stopped with status %d' % status)




