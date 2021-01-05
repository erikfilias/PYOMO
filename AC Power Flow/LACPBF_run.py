# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 17:34:21 2018

@author: Erik
"""

#from pyomo.environ import *
#from pyomo.opt import SolverFactory
from pyomo.environ import SolverFactory
from LACPBF import model
#from pyomo.core import Param
# ... more constraints
model.pprint()

instance = model.create('sistema14nos.dat')
# instance.pprint()

opt= SolverFactory("gurobi")

results = opt.solve(instance)

results.write()

print("------------------------------------------------------------------------------------------------------------")
print("-------------------------------SUMMARY----------------------------------------------------------------------")
print("------------------------------------------------------------------------------------------------------------")
print("BUS RESULTS")
print("------------------------------------------------------------------------------------------------------------")
print("   Bus      V[pu]  Theta[Degree]     Pg[MW]     Qg[MVAr]       Pd[MW]     Qd[MVAr]      Gsh[MW]    Bsh[MVAr]")
print("------------------------------------------------------------------------------------------------------------")
for i in instance.BAR:
    print('%5d  %10.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f'
          %(i,(instance.bus_voltage_sqr[i].value)**0.5,instance.bus_angle[i].value*180/3.14159265359,
            instance.Sbase*instance.bus_Pg[i].value,instance.Sbase*instance.bus_Qg[i].value,
            instance.Sbase*instance.Pd_pu[i],instance.Sbase*instance.Qd_pu[i],
            instance.Sbase*instance.gshb[i]*instance.bus_voltage_sqr[i].value,
            instance.Sbase*instance.bshb[i]*instance.bus_voltage_sqr[i].value))


a=0
b=0
c=0
d=0
e=0
f=0
for i in instance.BAR:
    a = a + instance.bus_Pg[i].value
    b = b + instance.bus_Qg[i].value
    c = c + instance.Pd[i]
    d = d + instance.Qd[i]
    e = e + instance.gshb[i]*instance.bus_voltage_sqr[i].value
    f = f + instance.bshb[i]*instance.bus_voltage_sqr[i].value
print('TOTAL %37.4f %12.4f %12.4f %12.4f %12.4f %12.4f'
      %(instance.Sbase*a,instance.Sbase*b,c,d,instance.Sbase*e,instance.Sbase*f))    

print("------------------------------------------------------------------------------------------------------------")
print("BRANCH RESULTS")
print("------------------------------------------------------------------------------------------------------------")
print("    i     j        Pij[MW]         Pji[MW]       Qij[MVAr]       Qji[MVAr]         Pls[MW]       Qls[MVAr]")
print("------------------------------------------------------------------------------------------------------------")
for i,j in instance.RAM:
	print('%5d %5d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f'
       %(i,j,
         instance.Sbase*(instance.branch_P[i,j].value+instance.r[i,j]*instance.branch_Iij_sqr[i,j].value),
         -instance.Sbase*instance.branch_P[i,j].value,
         instance.Sbase*(instance.branch_Q[i,j].value+instance.x[i,j]*instance.branch_Iij_sqr[i,j].value-instance.bsh_half[i,j]*instance.bus_voltage_sqr[i].value),
         -instance.Sbase*(instance.branch_Q[i,j].value+instance.bsh_half[i,j]*instance.bus_voltage_sqr[j].value),
         instance.Sbase*instance.r[i,j]*instance.branch_Iij_sqr[i,j].value,
         instance.Sbase*(instance.x[i,j]*instance.branch_Iij_sqr[i,j].value-instance.bsh_half[i,j]*instance.bus_voltage_sqr[i].value-instance.bsh_half[i,j]*instance.bus_voltage_sqr[j].value)))


a=0
b=0
for i,j in instance.RAM:
    a = a + instance.r[i,j]*instance.branch_Iij_sqr[i,j].value
    b = b + instance.x[i,j]*instance.branch_Iij_sqr[i,j].value-instance.bsh_half[i,j]*instance.bus_voltage_sqr[i].value-instance.bsh_half[i,j]*instance.bus_voltage_sqr[j].value
print("------------------------------------------------------------------------------------------------------------")
print('TOTAL %85.4f %15.4f'
      %(instance.Sbase*a,instance.Sbase*b))   
print("------------------------------------------------------------------------------------------------------------")
# for i,j in instance.RAM:
#     print('%15.8f %15.8f %15.8f %15.8f'
#           %(instance.z2[i,j],instance.r[i,j],instance.x[i,j],instance.bsh_half[i,j]))
# print("------------------------------------------------------------------------------------------------------------")
