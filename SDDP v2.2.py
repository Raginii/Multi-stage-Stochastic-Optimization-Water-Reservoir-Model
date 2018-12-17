# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 15:29:53 2018

@author: visakh
"""

from gurobipy import *
import matplotlib.pyplot as plt
import pandas as pd
import project_data as data
import random
import math

t = data.t
w = data.w
ft = data.ft 
h = data.h
wdemand = data.demandnew
wexchange = data.exchangeNew
wlist = data.wlist
capacity = data.wcapacity
flist = data.flist
wprices = data.wprices
winflow = data.winflow
swinflow = data.swinflow
resmax = data.resmax
resmin = data.resmin
s = data.slist
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
gencost = data.wgencost ## 
capacityKeys = capacity.keys()

w = w[:52]
s = s[:12]
t  = w
resmax = 106.2e6 #data.resmax
resmin = 10e6 #data.resmin
#wdemand = {a:0 for a in wlist}
#wexchange = {a:0 for a in wlist}
#wcapacity = {(a,f):0 for a in wlist for f in ft}
constrInflow = 0
constrDemand = 0
constrCapacity = {f:0 for f in ft}

#for (x,y) in wt:
#    wdemand[x] += demand[y]
#    wexchange[x] += exchange[y] 
    

##########################################################
#Creating Model
m = Model('SDDP')
m.params.logtoconsole=0

#Creating variables
gap = m.addVar(name = 'gap') ##Gap generation in T
res = m.addVar(name = 'res') ##reservoir energy level at end in T
spill = m.addVar(name = 'spill') ##Spillafe in hour t
x = m.addVars(ft,name = 'x') ##Generated power in hour t with ft
slackup = m.addVar(name = 'slackup') ##Slack for the upper level
slacklo = m.addVar(name = 'slacklo') ##Slack for the lower level
alpha = m.addVar(name = 'alpha')               

#update model
m.update()

#Creating constraints
hyd = m.addConstr(((res + x['Hydro'] + spill == constrInflow)) , "Hydraulic continuity" ) #inflow = winflow + resprev
dem = m.addConstr((quicksum(x[f] for f in ft)  + gap >= constrDemand) , "Demand Satisfaction" )
m.addConstr((-res + slackup >= -resmax) , "Maximum reservoir level" )
m.addConstr((res + slacklo >= resmin) , "Minimum reservoir level" )
cap = m.addConstrs(((x[f] <= constrCapacity[f]) for f in ft) , "Capacity limit" )

m.update()

converge = 1
iteration = 0
resPrevious={}
Constr = {i:[] for i in w}
lb_iter = {}
ub_iter = {}
sigma = {}
trialList={}

#clearing flow
while converge == 1 and iteration < 200:
    #Initializing iteration
    iteration = iteration + 1
    trialSet = random.sample(s,2)
#    print(trialSet)
    ub_iter_set = {}     
    #Forward propoation
    resPrev = resmin #Setting first week
    for week in w:
        if week == w[0]:
            newCuts = {}
            CutCount = 0
            hyd.RHS = winflow[week]*5+resPrev
            m.update()
            dem.RHS = wdemand[week] - wexchange[week]
            m.update()
            for f in ft:
                if (week,f) in capacityKeys:
                    cap[f].RHS = capacity[week,f]
                else:
                    cap[f].RHS = 0
            for i in Constr[week]: ##Adding benders cuts
                newCuts[CutCount]=m.addConstr(i[0],i[1],i[2])
                m.update()
                CutCount+=1
            m.setObjective((quicksum(gencost[week,f]*x[f] for f in ft)+ 1000*gap + 10e6*(slackup+slacklo)+ alpha),GRB.MINIMIZE)
            m.update()
            m.reset()
            m.optimize()
            for trial in trialSet:
                resPrevious[w[w.index(week)+1],trial] = res.X
            objVal = m.objval
            lb_iter[iteration]=objVal
            ub_iter_stageOne = (objVal - alpha.X)
#            print(iteration,week,ub_iter_stageOne,objVal,alpha.X)
            for i in range(CutCount): ##removing benders cuts
                m.remove(newCuts[i]) 
                m.update()
        else:
            for trial in trialSet:
                if (iteration,trial) not in ub_iter_set:
#                    print(iteration)
                    ub_iter_set[iteration,trial]=ub_iter_stageOne
                newCuts = {}
                CutCount = 0
                hyd.RHS = swinflow[trial,week]*2+3*winflow[week]+resPrevious[week,trial]
                m.update()
                dem.RHS = wdemand[week] - wexchange[week]
                m.update()
                for f in ft:
                    if (week,f) in capacityKeys:
                        cap[f].RHS = capacity[week,f]
                    else:
                        cap[f].RHS = 0
                for i in Constr[week]: ##Adding benders cuts
                    newCuts[CutCount]=m.addConstr(i[0],i[1],i[2])
                    m.update()
                    CutCount+=1
                m.setObjective((quicksum(gencost[week,f]*x[f] for f in ft)+ 1000*gap + 10e6*(slackup+slacklo)+ alpha),GRB.MINIMIZE)
                m.update()
                m.reset()
                m.optimize()
#                m.write("Check_%s_%s.lp"%(iteration,week))
                if week!=w[-1]:             
                    resPrevious[w[w.index(week)+1],trial] = res.X
                objVal = m.objval
#                print("Before",ub_iter_set[iteration,trial])
                if week==w[-1]:
                    ub_iter_set[iteration,trial]+= objVal
                else:
                    ub_iter_set[iteration,trial]+= (objVal - alpha.X)
#                print("After",ub_iter_set[iteration,trial])
#                print(iteration,week,trial,ub_iter_set[iteration,trial],objVal,alpha.X)
#                print(iteration,trial,week,(objVal - alpha.X))
                for i in range(CutCount): ##removing benders cuts
                    m.remove(newCuts[i]) 
                    m.update()
#    print(iteration,ub_iter_set)
#    print(ub_iter_set)
    ub_iter[iteration]=sum(ub_iter_set[iteration,trial] for trial in trialSet)/len(trialSet)
    #Backward Propogation
    for trial in trialSet:
        for week in w[::-1]:
            v=0
            h=0
            for scen in s:
                CutCount=0
                m.update()
                m.reset()
                #m.remove(hyd)
                for i in Constr[week]: ##Adding benders cuts
                    newCuts[CutCount]=m.addConstr(i[0],i[1],i[2])
                    m.update()
                    CutCount+=1
                if week == w[0]:
                    #hyd = m.addConstr(((res + x['Hydro'] + spill <= winflow[week]+resPrev)) , "Hydraulic continuity" )
                    hyd.RHS = winflow[week]*5+resPrev
                    m.update()
                else:
                    #hyd = m.addConstr(((res + x['Hydro'] + spill <= swinflow[scen,week]+resPrev)) , "Hydraulic continuity" )
                    hyd.RHS = swinflow[scen,week]*2+3*winflow[week]+resPrevious[week,trial]
                    m.update()
                dem.RHS = wdemand[week] - wexchange[week]
                m.update()
                for f in ft:
                    if (week,f) in capacityKeys:
                        cap[f].RHS = capacity[week,f]
                        m.update()
                    else:
                        cap[f].RHS = 0
                        m.update()
                m.setObjective((quicksum(gencost[week,f]*x[f] for f in ft)+ 1000*gap + 10e6*(slackup+slacklo)+ alpha),GRB.MINIMIZE)                
                m.update()
                m.reset()
                m.optimize()
#                m.write("back_%s_%s_%s.lp"%(iteration,week,scen))
                #m.write("backDetail_%s_%s_%s.mps"%(iteration,week,scen))
                #print("Duals",hyd.Pi,dem.Pi)
                lamdaVal = hyd.Pi
                #print(lamdaVal,hyd.Pi,res.X)
                #print("Backward",iteration,week,scen,hyd.Pi,res.X,m.objVal)
#                print(iteration,week,scen,"Hydro",x['Hydro'].X)
                objVal = m.objval
                m.reset()
                for i in range(CutCount): ##removing benders cuts
                    m.remove(newCuts[i]) 
                    m.update()
        #                print(iteration,week,scen,lamdaVal)
                if week==w[0]:
                    v+=(1/len(s))*(objVal - (lamdaVal*resPrev))
                else:
                    v+=(1/len(s))*(objVal - (lamdaVal*resPrevious[week,trial]))
                h+=(1/len(s))*lamdaVal
            #Evaluating Benders                
            if week != w[0]:
                ConstrNew = m.addConstr(alpha >= v + (h*res))
                m.update()
                lhs = m.getRow(ConstrNew)
                sense = ConstrNew.Sense
                rhs = ConstrNew.RHS
                Constr[w[w.index(week)-1]].append((lhs,sense,rhs))
                m.remove(ConstrNew)
                m.update()
    sigma[iteration]= math.sqrt((sum(((ub_iter[iteration]-ub_iter_set[iteration,trial])**2) for trial in trialSet))/(len(trialSet)**2))
    if (lb_iter[iteration]>=(ub_iter[iteration]-(2*sigma[iteration]))) and (lb_iter[iteration]<=(ub_iter[iteration]+(2*sigma[iteration]))):
        print("Convergence reached at iteration",iteration)
        print(lb_iter[iteration],ub_iter[iteration])
        converge=0

lb = lb_iter[iteration]
#for i in range(1,iteration+1):
#    print(i,"lower bound:",lb_iter[i],"upper bound: ",ub_iter[i],"\n")
print("lower bound current",lb)    
#Upper bound simulation
k=9
trialSet= [i for i in range(k)]
while converge == 1 and iteration < 2000:
    #Initializing iteration
#    print(iteration)
    iteration = iteration + 1
    lb_iter[iteration]=lb
    ub_iter_set = {}     
    #Forward propoation
    resPrev = resmin #Setting first week
    for week in w:
        if week == w[0]:
            newCuts = {}
            CutCount = 0
            hyd.RHS = winflow[week]*5+resPrev
            m.update()
            dem.RHS = wdemand[week] - wexchange[week]
            m.update()
            for f in ft:
                if (week,f) in capacityKeys:
                    cap[f].RHS = capacity[week,f]
                else:
                    cap[f].RHS = 0
            for i in Constr[week]: ##Adding benders cuts
                newCuts[CutCount]=m.addConstr(i[0],i[1],i[2])
                m.update()
                CutCount+=1
            m.setObjective((quicksum(gencost[week,f]*x[f] for f in ft)+ 1000*gap + 10e6*(slackup+slacklo)+ alpha),GRB.MINIMIZE)
            m.update()
            m.reset()
            m.optimize()
            for trial in trialSet:
                resPrevious[w[w.index(week)+1],trial] = res.X
            objVal = m.objval
            ub_iter_stageOne = (objVal - alpha.X)
            for i in range(CutCount): ##removing benders cuts
                m.remove(newCuts[i]) 
                m.update()
        else:
            for trial in trialSet:
                trialList[trial,week]=random.choice(s)
                if (iteration,trial) not in ub_iter_set:
                    ub_iter_set[iteration,trial]=ub_iter_stageOne
                newCuts = {}
                CutCount = 0
                hyd.RHS = swinflow[trialList[trial,week],week]*2+3*winflow[week]+resPrevious[week,trial]
                m.update()
                dem.RHS = wdemand[week] - wexchange[week]
                m.update()
                for f in ft:
                    if (week,f) in capacityKeys:
                        cap[f].RHS = capacity[week,f]
                    else:
                        cap[f].RHS = 0
                for i in Constr[week]: ##Adding benders cuts
                    newCuts[CutCount]=m.addConstr(i[0],i[1],i[2])
                    m.update()
                    CutCount+=1
                m.setObjective((quicksum(gencost[week,f]*x[f] for f in ft)+ 1000*gap + 10e6*(slackup+slacklo)+ alpha),GRB.MINIMIZE)
                m.update()
                m.reset()
                m.optimize()
                if week!=w[-1]:             
                    resPrevious[w[w.index(week)+1],trial] = res.X
                objVal = m.objval
                alphaValue = alpha.X
                ub_iter_set[iteration,trial]+= (objVal - alphaValue)
                for i in range(CutCount): ##removing benders cuts
                    m.remove(newCuts[i]) 
                    m.update()
    ub_iter[iteration]=sum(ub_iter_set[iteration,trial] for trial in trialSet)/len(trialSet)
    sigma[iteration]= math.sqrt((sum(((ub_iter[iteration]-ub_iter_set[iteration,trial])**2) for trial in trialSet))/(len(trialSet)**2))
    if (lb_iter[iteration]>=(ub_iter[iteration]-(1.96*sigma[iteration]))) and (lb_iter[iteration]<=(ub_iter[iteration]+(1.96*sigma[iteration]))):
        print("Convergence reached at iteration",iteration)
        print(lb_iter[iteration],ub_iter[iteration])
        converge=0            
if converge!=0:
    print("End Bounds",lb,ub_iter[iteration])      

plt.plot(lb_iter.keys(),lb_iter.values())
plt.plot(ub_iter.keys(),ub_iter.values())                            
                
    

