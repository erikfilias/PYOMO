#from pyomo.environ import *
#from pyomo.opt import SolverFactory
from pyomo.environ import SolverFactory
import pyomo.environ as pe
# import ACPFR
# from PYOMO import ACPFR

# include ACPFR

# include sistema14nos.dat
# load "sistema14nos"
#from pyomo.core import Param
# ... more constraints
model.pprint()

instance = model.create_instance('sistema14nos.dat')
instance.pprint()

opt= SolverFactory('ipopt')

results = opt.solve(instance)

results.write()

print("------------------------------------------------------------------------------------------------------------")
print("-------------------------------SUMMARY----------------------------------------------------------------------")
print("------------------------------------------------------------------------------------------------------------")
print("\nBUS RESULTS\n")
print("------------------------------------------------------------------------------------------------------------")
print("   Bus      V[pu]  Theta[Degree]     Pg[MW]     Qg[MVAr]       Pd[MW]     Qd[MVAr]      Gsh[MW]    Bsh[MVAr]\n")
#for i in range(len(model.bus_angle)): print("Bus_Angle "+str(i)+":" + str(model.bus_angle[i].value*180/3.14159265359))
for i in instance.BAR:
#    print('%5d  %10.4f %12.4f %12.4f %12.4f '%(i,(instance.bus_voltage_sqr[i].value)**0.5,instance.bus_angle[i].value*180/3.14159265359,100*instance.bus_Pg[i].value,100*instance.bus_Qg[i].value))
    a = (instance.bus_e[i].value**2 + instance.bus_f[i].value**2)
    b = 180/3.14159265359*pe.atan((instance.bus_f[i].value)/instance.bus_e[i].value)
    print('%5d  %10.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f'
          %(i,a**0.5,b,
            instance.Sbase*instance.bus_Pg[i].value,instance.Sbase*instance.bus_Qg[i].value,
            instance.Sbase*instance.Pd_pu[i],instance.Sbase*instance.Qd_pu[i],
            instance.Sbase*instance.gshb[i]*a,
            instance.Sbase*instance.bshb[i]*a))


a=0
b=0
c=0
d=0
e=0
f=0
for i in instance.BAR:
    a1 = (instance.bus_e[i].value**2 + instance.bus_f[i].value**2)
    a = a + instance.bus_Pg[i].value
    b = b + instance.bus_Qg[i].value
    c = c + instance.Pd[i]
    d = d + instance.Qd[i]
    e = e + instance.gshb[i]*a1
    f = f + instance.bshb[i]*a1
print('TOTAL %37.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n'
      %(instance.Sbase*a,instance.Sbase*b,c,d,instance.Sbase*e,instance.Sbase*f))
#
print("------------------------------------------------------------------------------------------------------------")
print("\nBRANCH RESULTS\n")
print("------------------------------------------------------------------------------------------------------------")
print("    i     j        Pij[MW]         Pji[MW]       Qij[MVAr]       Qji[MVAr]         Pls[MW]       Qls[MVAr]\n")
for i,j in instance.RAM:
	print('%5d %5d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f'
       %(i,j,
         instance.Sbase*(instance.branch_Pde[i,j].value),
         instance.Sbase*instance.branch_Ppara[i,j].value,
         instance.Sbase*(instance.branch_Qde[i,j].value),
         instance.Sbase*(instance.branch_Qpara[i,j].value),
         ((instance.Sbase*(instance.branch_Pde[i,j].value+instance.branch_Ppara[i,j].value))**2)**0.5,
         ((instance.Sbase*(instance.branch_Qde[i,j].value+instance.branch_Qpara[i,j].value))**2)**0.5))



a=0
b=0
for i,j in instance.RAM:
    a = a + instance.branch_Ppara[i,j].value+instance.branch_Pde[i,j].value
    b = b + instance.branch_Qpara[i,j].value+instance.branch_Qde[i,j].value
print('TOTAL %85.4f %15.4f\n'
      %(instance.Sbase*a,instance.Sbase*b))
