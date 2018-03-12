import numpy as np
from math import pi
import matplotlib.pyplot as plt
import numpy as np
#EVERYTHING IN SI UNITS#########



#variables
S = 30.
e = 0.8
C_Lalpha = 5.084
A = 8.438664
cd0t = 0.04
cd0new = 0.0228462731283
alpha0 = -1.109998

payload = 82+95+67+84+60+73+78+74+90+89
fuelweight = 2560*0.45359237
BEW = 9165*0.45359237
rampweight = BEW + payload + fuelweight

Press0 = 101325.
Temp0 = 288.15
Rho0 = 1.225
Lambda = -0.0065
R = 287.05
g = 9.80665
gamma = 1.4

def weight(Wf):
    return g*(rampweight-Wf)

def pressure(Alt):    # Static Pressure based on pressure altitude
    return Press0 * (1 + Lambda * Alt / Temp0)**(-g/(Lambda*R))   # In [Pa]

def Mach(Pressure,Vc): # Mach Number based on calibrated airspeed
    return np.sqrt(2./(gamma-1.)*((1.+(Press0/Pressure)*((1. + (gamma - 1.)/ (2.*gamma) * Rho0 / Press0 * Vc ** 2)**(gamma / (gamma-1.))-1.))**((gamma-1.)/gamma) - 1.)) # Dimensionless

def temp(TAT,Mach): # Static Temperature
    return TAT / (1. + (gamma-1.)/2. * Mach ** 2) # In [K], modified for ram rise

def difftempISA(temp,alt):
    return temp-(Temp0 + Lambda*alt)

def speedofsound(temp):
    return np.sqrt(gamma * R * temp)  # In [m/s]

def TAS(Mach,a): # True Airspeed (corrected for ram rise)
    return Mach * a # in [m/s] also

def density(pressure,temp): # Air Density
    return pressure / (temp * R) # in [kg/m^3]

def EAS(TAS,density): # Equivalent Airspeed (based on calibrated airspeed)
    return TAS * np.sqrt(density / Rho0)   # in [m/s]

def cl(weight,density,TAS):
    return weight / (0.5 * density * TAS**2 * S)

def cd(thrust,density,TAS): # Drag Coefficient
    return thrust / (0.5 * density * TAS** 2 * S)

def cd0(cd,cl): # Parasitic Drag
    return cd - cl **2 / (pi * e * A)

hft = [11990,12010,12000,12000,11970,11980]
alt = 6*[0]
for i in range(len(hft)):
    alt[i] = hft[i] * 0.3048

IASkts = [245,224,191,160,128,112]
IAS = 6*[0]
for i in range(len(hft)):
    IAS[i] = IASkts[i] *  0.514444444

AoA = [1.5,1.9,3.2,4.9,8.2,11]

cltheo = 6 * [0]
for i in range(len(AoA)):
    cltheo[i] = C_Lalpha*np.radians((AoA[i] - alpha0))

TATC = [-1.4,-3.4,-5.6,-7.2,-9.2,-9.2]
TAT = 6*[0]
for i in range(len(TATC)):
    TAT[i] = TATC[i] + 273.15

thrust = [(3460.78	+ 3610.77),(2867.75	+3038.24),(2187.26	+2300.58),(1973.54	+2165.12),(1678.75	+1871.05),(2028.88	+2195.41)]

Wflbs = [498,538,592,620,640,672]
Wf = 6*[0]
for i in range(len(Wflbs)):
    Wf[i] = Wflbs[i] * 0.45359237

result = [[],[],[]]

for i in range(6):
    pressure2 = pressure(alt[i])
    mach2 = Mach(pressure2,IAS[i])
    temp2 = temp(TAT[i],mach2)
    density2 = density(pressure2,temp2)
    a2 = speedofsound(temp2)
    TAS2 = TAS(mach2,a2)
    EAS2 = EAS(TAS2,density2)
    weight2 = weight(Wf[i])
    Cl2 = cl(weight2,density2,TAS2)
    Cd2 = cd(thrust[i],density2,TAS2)
    Cd02 = cd0(Cd2,Cl2)
    #print('pressure: {}, mach: {}, temp: {}, density: {}, speed of sound: {}, TAS: {}, EAS: {}, Weight: {}, Cl: {}, CD: {}, Cd0: {}'.format(pressure2,mach2,temp2,density2,a2,TAS2,EAS2,weight2,Cl2,Cd2,Cd02))
    result[0].append(AoA[i])
    result[1].append(Cl2)
    result[2].append(Cd2)

cd0new = 0
for i in range(6):
    cd0new += cd0(result[2][i],result[1][i])
cd0new = cd0new/6.

cdtheo = 6*[0]
for i in range(len(cdtheo)):
    cdtheo[i] = cd0t + result[1][i]**2/(pi*e*A)
cdtheo2 = 6*[0]
for i in range(len(cdtheo2)):
    cdtheo2[i] = cd0new + result[1][i]**2/(pi*e*A)

#CL-ALPHA PLOT
x,y = result[0],result[1]
x2,y2 = result[0],cltheo
plt.scatter(x,y,label='practical')
plt.plot(x2,y2,label='theoretical')
# calc the trendline
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
plt.grid()
plt.ylabel('CL [-]')
plt.xlabel('Angle of attack [degrees]')
plt.axis([0, 15, 0, 1.25])
plt.legend()
plt.show()

#CL-ALPHA PLOT
x,y = result[0],result[2]
x2,y2 = result[0],cdtheo
x3,y3 = result[0],cdtheo2
plt.scatter(x,y,label='practical')
plt.plot(x2,y2,label='theoretical')
plt.plot(x3,y3,label='theoretical (with new CD0)')
# calc the trendline
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
plt.grid()
plt.ylabel('CD [-]')
plt.xlabel('Angle of attack [degrees]')
plt.axis([0, 15, 0, 0.1])
plt.legend()
plt.show()

#CL-CD PLOT
x,y = result[2],result[1]
x2,y2 = cdtheo,cltheo
x3,y3 = cdtheo2,cltheo
plt.scatter(x,y,label='practical')
plt.plot(x2,y2,label='theoretical')
plt.plot(x3,y3,label='theoretical (with new CD0)')
# calc the trendline
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
plt.grid()
plt.ylabel('CL [-]')
plt.xlabel('CD [-]')
plt.axis([0, 0.1, 0, 1.2])
plt.legend()
plt.show()

clsquare = 6*[0]
for i in range(len(clsquare)):
    clsquare[i] = result[1][i]**2

cltsquare = 6*[0]
for i in range(len(cltsquare)):
    cltsquare[i] = cltheo[i]**2

#CD-CL^2 PLOT
x,y = clsquare,result[2]
x2,y2 = cltsquare, cdtheo
x3,y3 = cltsquare, cdtheo2
plt.scatter(x,y,label='practical')
plt.plot(x2,y2,label='theoretical')
plt.plot(x3,y3,label='theoretical (with new CD0)')
# calc the trendline
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
plt.grid()
plt.ylabel('CD [-]')
plt.xlabel('CL^2 [-]')
plt.axis([0, 1.2, 0, 0.1])
plt.legend()
plt.show()