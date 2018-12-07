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
capacity = data.wcapacity
gencost = data.wgencost ## 
capacityKeys = capacity.keys()

w = w[:8]
s = s[:3]
t  = w
resmax = 106.2e6 #data.resmax
resmin = 10e6 #data.resmin
wdemand = {a:0 for a in wlist}
wexchange = {a:0 for a in wlist}

for (x,y) in wt:
    wdemand[x] += demand[y]
    wexchange[x] += exchange[y] 
   

print(winflow)
##########################################################
#Creating Model
m = Model('SDDP_Hydro')

#Creating variables
gap = m.addVars(t,name = 'gap') ##Gap generation in T
res = m.addVars(t,name = 'res') ##reservoir energy level at end in T
spill = m.addVars(t,name = 'spill') ##Spillafe in hour t
x = m.addVars(t,ft,name = 'x') ##Generated power in hour t with ft
slackup = m.addVars(t,name = 'slackup') ##Slack for the upper level
slacklo = m.addVars(t,name = 'slacklo') ##Slack for the lower level
alpha = m.addVars(w,name = 'alpha')               
#update model
m.update()

#Creating constraints
hyd = m.addConstrs(((res[i] - res[t[t.index(i)-1]] + x[i,'Hydro'] + spill[i] == winflow[i]) for i in t if i!= 'w1') , "Hydraulic continuity" )
m.addConstrs(((quicksum(x[i,f] for f in ft if (i,f) in capacityKeys)  + gap[i] >= wdemand[i] - wexchange[i]) for i in t) , "Demand Satisfaction" )
m.addConstrs(((-res[i] + slackup[i] >= -resmax) for i in t) , "Maximum reservoir level" )
m.addConstrs(((res[i] + slacklo[i] >= resmin) for i in t) , "Minimum reservoir level" )
m.addConstrs(((-x[i,f] >= -capacity[i,f]) for i in t for f in ft if (i,f) in capacityKeys) , "Minimum reservoir level" )
m.addConstr(((res['w1'] - resmin + x['w1','Hydro'] + spill['w1'] == winflow['w1']) ) , "Hydraulic continuity" )
#Setting objective
m.setObjective(
        (
        quicksum(gencost[i,f]*x[i,f] for i in t for f in ft if (i,f) in capacityKeys)
        + 1000*quicksum(gap[i] for i in t)
        + 10e6*quicksum(slackup[i]+slacklo[i] for i in t)+ quicksum(alpha[i] for i in w)
        )
        ,GRB.MINIMIZE)

m.update()


converge = 1
resVal = {}
objVal = {}
iteration = 0

while converge == 1 and iteration < 100:
	iteration + =1
	Cons = {}
	for weeks in w:
		# Optimize and get the reservoir level values per week 
		# (Used deterministic inflow)
		m.optimize()
		resVal[weeks] = res[weeks].X
		objVal[weeks] = m.objVal

		# Fix res level using the previous week
		Cons[weeks] = m.addConstr(res[weeks] == resVal[weeks],"ForwardPassFix")
	# for i in w:
	# 	m.remove(Cons[i])

	# Set the lower bound to week 1's objective value.
	lb = objVal['w1']
	

	### Backward Pass
	for weeks in w[::-1]:
		for scen in s:
			m.optimize()
			for i in hyd:
					i.RHS = 








