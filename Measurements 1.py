# Import necessary modules
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import numpy as np

# EVERYTHING IN SI UNITS#########

# Variables
S = 30.                 # Wing area
e = 0.8                 # Given Oswald efficiency factor
C_Lalpha = 5.084        # 4.4284 is the measured clalpha.
A = 8.438664            # Aspect ratio
cd0t = 0.04             # Theoretical zero-lift drag coefficient
cd0new = 0.0228462731283 # Measured zero-lift drag coefficient
alpha0 = -1.109998      # Measured alpha0
chord = 2.0569

payload = 82+95+67+84+60+73+78+74+90+89 # Total payload mass (passengers)
fuelweight = 2560*0.45359237            # Fuel mass at the start
BEW = 9165*0.45359237                   # Basic empty mass
rampweight = BEW + payload + fuelweight # Total mass

#ISA Properties
Press0 = 101325.
Temp0 = 288.15
Rho0 = 1.225
Lambda = -0.0065
R = 287.05
g = 9.80665
gamma = 1.4
beta = 1.458*10**-6  # constant to calculate dynamic viscosity of air
Sair= 110.4   # constant to calculate dynamic viscosity of air

def weight(Wf):     # Calculate the weight at a given fuelweight used
    return g*(rampweight-Wf)

def pressure(Alt):    # Static Pressure based on pressure altitude
    return Press0 * (1 + Lambda * Alt / Temp0)**(-g/(Lambda*R))   # In [Pa]

def Mach(Pressure,Vc): # Mach Number based on calibrated airspeed
    return np.sqrt(2./(gamma-1.)*((1.+(Press0/Pressure)*((1. + (gamma - 1.)/ (2.*gamma) * Rho0 / Press0 * Vc ** 2)**(gamma / (gamma-1.))-1.))**((gamma-1.)/gamma) - 1.)) # Dimensionless

def temp(TAT,Mach): # Static Temperature
    return TAT / (1. + (gamma-1.)/2. * Mach ** 2) # In [K], modified for ram rise

def difftempISA(temp,alt): # Delta T
    return temp-(Temp0 + Lambda*alt)

def speedofsound(temp): # Speed of sound
    return np.sqrt(gamma * R * temp)  # In [m/s]

def TAS(Mach,a): # True Airspeed (corrected for ram rise)
    return Mach * a # in [m/s] also

def density(pressure,temp): # Air Density
    return pressure / (temp * R) # in [kg/m^3]

def EAS(TAS,density): # Equivalent Airspeed (based on calibrated airspeed)
    return TAS * np.sqrt(density / Rho0)   # in [m/s]

def cl(weight,density,TAS): # Lift coefficient
    return weight / (0.5 * density * TAS**2 * S)

def cd(thrust,density,TAS): # Drag Coefficient
    return thrust / (0.5 * density * TAS** 2 * S)

def cd0(cd,cl): # Parasitic Drag
    return cd - cl **2 / (pi * e * A)


def viscosity(temp):  # dynamic viscosity of air given temperature
    return beta * temp ** (3 / 2) / (temp + Sair)


def reynolds(density, velocity, temperature):  # reynolds number
    return density * velocity * chord / viscosity(temperature)

hft = [11990,12010,12000,12000,11970,11980] # Different altitudes for the measurements, corrected to SI units
alt = 6*[0]
for i in range(len(hft)):
    alt[i] = hft[i] * 0.3048

IASkts = [245,224,191,160,128,112] # Different velocities for the measurements, corrected to SI units, converting to calibrated airspeed
IAS = 6*[0]
for i in range(len(hft)):
    IAS[i] = (IASkts[i]-2.) *  0.514444444

AoA = [1.5,1.9,3.2,4.9,8.2,11] # Different angles of attack for the measurements

cltheo = 6 * [0]                # Theoretical lift coefficient
for i in range(len(AoA)):
    cltheo[i] = C_Lalpha*np.radians((AoA[i] - alpha0))

TATC = [-1.4,-3.4,-5.6,-7.2,-9.2,-9.2]      # Different temperatures for the measurements, corrected to SI units
TAT = 6*[0]
for i in range(len(TATC)):
    TAT[i] = TATC[i] + 273.15

# Thrust following from the thruxt.exe file
thrust = [(3460.78	+ 3610.77),(2867.75	+3038.24),(2187.26	+2300.58),(1973.54	+2165.12),(1678.75	+1871.05),(2028.88	+2195.41)]

Wflbs = [498,538,592,620,640,672]  # # Different weights of fuel used for the measurements, corrected to SI units
Wf = 6*[0]
for i in range(len(Wflbs)):
    Wf[i] = Wflbs[i] * 0.45359237

result = [[],[],[]]
reynoldslist = 6*[0]
machlist = 6*[0]

for i in range(6): # Calculating all parameters for all measurement points
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
    reynoldslist[i] = reynolds(density2,TAS2,temp2)
    machlist[i] = mach2

cd0new = 0 # Zero-lift drag coefficient
for i in range(6):
    cd0new += cd0(result[2][i],result[1][i])
cd0new = cd0new/6.

cdtheo = 6*[0] # Theoretical drag coefficient
for i in range(len(cdtheo)):
    cdtheo[i] = cd0t + result[1][i]**2/(pi*e*A)
cdtheo2 = 6*[0]
for i in range(len(cdtheo2)):
    cdtheo2[i] = cd0new + result[1][i]**2/(pi*e*A)



# Set plot size
plt.rcParams.update({'font.size': 15})


#CL-ALPHA PLOT
x,y = result[0],result[1]
x2,y2 = result[0],cltheo
plt.scatter(x,y,label='measurement points')
plt.plot(x2,y2,linestyle=':',label='theoretical')
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--",label='polynomial fit')
plt.grid()
plt.ylabel('CL [-]')
plt.xlabel('Angle of attack [degrees]')
plt.axis([0, 12, 0, 1.5])
plt.legend(loc=2)
plt.show()

#CL-ALPHA PLOT
x,y = result[0],result[2]
x2,y2 = result[0],cdtheo
x3,y3 = result[0],cdtheo2
plt.scatter(x,y,label='measurement points')
plt.plot(x2,y2,linestyle=':',label='theoretical')
plt.plot(x3,y3,label='theoretical (new CD0)')
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
plt.plot(x,p(x),"r--",label='polynomial fit')
plt.grid()
plt.ylabel('CD [-]')
plt.xlabel('Angle of attack [degrees]')
plt.axis([0, 12, 0, 0.1])
plt.legend(loc=2)
plt.show()

Cl2nd = [0.4635,0.5267,0.5955,0.6892,0.3991,0.3595,0.3223]
Cd2nd = [0.031134,0.035899,0.04164,0.04916,0.02658,0.02372,0.021222]

#CL-CD PLOT
x,y = result[2],result[1]
x2,y2 = cdtheo,cltheo
x3,y3 = cdtheo2,cltheo
plt.scatter(x,y,label='measurement points')
plt.scatter(Cd2nd,Cl2nd,linestyle=':',label='2nd measurement')
#plt.plot(x2,y2,linestyle=':',label='theoretical')
plt.plot(x3,y3,label='theoretical (new CD0)')
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
plt.plot(x,p(x),"r--",label='polynomial fit')
plt.grid()
plt.ylabel('CL [-]')
plt.xlabel('CD [-]')
plt.axis([0, 0.1, 0, 1.4])
plt.legend(loc=2)
plt.show()

clsquare = 6*[0]  # Calculating the lift coefficient squared for the measurements
for i in range(len(clsquare)):
    clsquare[i] = result[1][i]**2

cltsquare = 6*[0] # Calculating the lift coefficient squared for the theoretical plot
for i in range(len(cltsquare)):
    cltsquare[i] = cltheo[i]**2

#CD-CL^2 PLOT
x,y = clsquare,result[2]
x2,y2 = cltsquare, cdtheo
x3,y3 = cltsquare, cdtheo2
plt.scatter(x,y,label='measurement points')
plt.plot(x2,y2,linestyle=':',label='theoretical')
plt.plot(x3,y3,label='theoretical (new CD0)')
# calc the trendline
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--",label='polynomial fit')
plt.grid()
plt.ylabel('CD [-]')
plt.xlabel('CL^2 [-]')
plt.axis([0, 1.2, 0, 0.1])
plt.legend(loc=2)
plt.show()