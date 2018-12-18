from gurobipy import *
import pandas as pd
import time
import project_data as data

t = data.t
w = data.w
ft = data.ft 
h = data.h
demand = data.demand
exchange = data.exchange
wlist = data.wlist
flist = data.flist
wprices = data.wprices
winflow = data.winflow
swinflow = data.swinflow
resmax = 106.2e6 #data.resmax
resmin = 10e6 #data.resmin

wt = data.wt

demandNew = {a:0 for a in wlist}
exchangeNew = {a:0 for a in wlist}

for (x,y) in wt:
    demandNew[x] += demand[y]
    exchangeNew[x] += exchange[y] 
s = data.slist # Contains the list pf scenarios

inflow = data.inflow  ## Not used 
capacity = data.wcapacity
gencost = data.wgencost  ## Weekly 
capacityKeys = capacity.keys()


n= [0] # Node List
timeN={}
ancestorN= {}
scenN = {}
j=0
root = [0]
p = {}
####
#w = w[:6]
w = w[:6]
s = s[:12]
p[0] = 1
timeN[0] = w[0]
for weeks in w:
    prev = root
    if weeks!= w[0]:
        root = []
        for nodes in prev:
            for scen in s:
                j+=1
                n.append(j)
                timeN[j] = weeks
                ancestorN[j] = nodes
                scenN[j] = scen
                p[j] = (1/len(s))*p[nodes]
                root.append(j)

Nlist  = n +[0]

##########################################################
#Creating Model
m = Model('NodeDam')
start = time.time()

#Creating variables
gap = m.addVars(n,name = 'gap') ##Gap generation in T
res = m.addVars(n,name = 'res') ##reservoir energy level at end in T
spill = m.addVars(n,name = 'spill') ##Spillafe in hour t
x = m.addVars(n,ft,name = 'x') ##Generated power in hour t with ft
slackup = m.addVars(n,name = 'slackup') ##Slack for the upper level
slacklo = m.addVars(n,name = 'slacklo') ##Slack for the lower level
                   
#update model
m.update()

#Creating constraints
m.addConstrs(((res[i] - res[ancestorN[i]] + x[i,'Hydro'] + spill[i] == swinflow[scenN[i],timeN[i]]) for i in  n if i!=0) , "Hydraulic continuity" )
m.addConstrs(((quicksum(x[i,f] for f in ft if (timeN[i],f) in capacityKeys)  + gap[i] >= demandNew[timeN[i]] - exchangeNew[timeN[i]]) for i in n) , "Demand Satisfaction" )
m.addConstrs(((-res[i] + slackup[i] >= -resmax) for i in n) , "Maximum reservoir level" )
m.addConstrs(((res[i] + slacklo[i] >= resmin) for i in n) , "Minimum reservoir level" )
m.addConstrs(((-x[i,f] >= -capacity[timeN[i],f]) for i in n for f in ft if (timeN[i],f) in capacityKeys) , "Minimum reservoir level" )
m.addConstr((res[0] - resmin + x[0,'Hydro'] + spill[0] == winflow[w[0]]*5) , "Hydraulic continuity For Node1" )

#Setting objective
m.setObjective(
        (
        quicksum(gencost[timeN[i],f]*x[i,f]*p[i] for i in n for f in ft if (timeN[i],f) in capacityKeys)
        + 1000*quicksum(gap[i]*p[i] for i in n)
        + 10e6*quicksum(slackup[i]*p[i]+slacklo[i]*p[i] for i in n)
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
totTime = time.time() - start
print("The node based model solved in : {:.2f}".format(totTime),"seconds")