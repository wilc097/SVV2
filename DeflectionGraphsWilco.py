import numpy as np
from math import pi
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as poly


#EVERYTHING IN SI UNITS#########

#variables
S = 30.
e = 0.8
C_Lalpha = 5.084
A = 8.438664
cd0 = 0.04
Ws = 60500. #standard aircraft mass

payload = 82+95+67+84+60+73+78+74+90+89
fuelweight = 2560*0.45359237
BEW = 9165*0.45359237
rampweight = BEW + payload + fuelweight

Press0 = 101325.
Temp0 = 288.15
Rho0 = 1.225
Lambda = -0.0065
R = 287.05
g = 9.81
gamma = 1.4

###
beta = 1.458*10**-6  # (constant to calculate dynamic viscosity of air)
Sair= 110.4   # (constant to calculate dynamic viscosity of air)
chord = 2.0569
Cmtc = -0.064 # dimensionless thrust moment arm
xcg1 = 0.472
xcg2 = 0.4244
deltaxcg= xcg2- xcg1
cmdtheory = -1.1642
cmalphatheory = -0.5626

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

def column(matrix, i): #returns column  #i of data from matrix
    return [row[i] for row in matrix]

def redEAS(eas,weight): #reduced equivalent airspeed - weight = weight of aircraft
    return eas * (Ws/weight)**(0.5)

def viscosity(temp): #dynamic viscosity of air given temperature
    return beta*temp**(3/2)/(temp+Sair)
    
def reynolds(density,velocity,temperature): #reynolds number
    return density*velocity*chord/viscosity(temperature) 

def redelevdef(elevatdefl,tcs,tc,Cmdelta): #measured elevator deflection, standardised thrust coeff, thrust coeff
        return elevatdefl - Cmtc*(tcs-tc)/Cmdelta

def redelvforce(fmeasured, weightfuelused): #reduced elevator force - weight of fuel used converted into weight of aircraft
    return fmeasured*Ws/weight(weightfuelused)

def Cmdelta(CN,changeXcg,changeelevdefl): #CN assumed to be CL
    return - CN*changeXcg/chord/changeelevdefl

def thrustcoeff(thrust,density,velocity,diameter):
    return thrust/(density*velocity**2*diameter**2)

### Second series of measurements data processing

# Inputting data for 2nd series
datadefl = np.genfromtxt('ElevatorDeflSI.txt', dtype=float, comments='#',  delimiter="&", skip_header=1)


#finding pressure
press = list()
for i in column(datadefl,0):
        press.append(pressure(i))
#print(press)

#finding mach
mach = list()
for i,x in zip(column(datadefl,1), press):
    mach.append(Mach(x,i))
#print(mach)    

#temperature
temper = list()
for i,x in zip(column(datadefl,9), mach):
    temper.append(temp(i,x))
#print(temper)    

#speed sound
ssound = list()
for i in temper:
    ssound.append(speedofsound(i))
#print(ssound)  

#density
dens = list()
for i,x in zip(press,temper):
    dens.append(density(i,x))
#print(dens)

#equivalent and true velocity
tas = list()
veq = list()
for i,x,y in zip(mach,ssound,dens):
    tas.append(TAS(i,x))
    veq.append(EAS(TAS(i,x),y))
#print(tas)
#print(veq)
    
#reduced eq velocity
redveq = list()
for i,x in zip(veq,column(datadefl,8)):
    redveq.append(redEAS(i,weight(x)))
#print(redveq)
sortedredveq = sorted(redveq, key=float) 

#temperature difference
tempdiff = list()
for i,x in zip(temper,column(datadefl,0)):
    tempdiff.append(difftempISA(i,x))
#print(tempdiff)

#Generating data for thrust.exe thrust and standardised thrust
dataforthrust, dataforthruststd = np.array([]),np.array([])
fuelflowstd = [0.048,0.048,0.048,0.048,0.048,0.048,0.048]
dataforthrust = np.dstack((column(datadefl,0),mach,tempdiff,column(datadefl,6),column(datadefl,7)))
dataforthruststd = np.dstack((column(datadefl,0),mach,tempdiff,fuelflowstd,fuelflowstd))
#print(dataforthruststd)



#getting  thrust data processed from thrust.exe
thrustLR2 = np.genfromtxt('thrust2.txt', dtype=float)
thrustLRstd2 = np.genfromtxt('thruststd2.txt', dtype=float)
thrust2, thruststd2= list(),list()
for i in range(7):
    t1=sum(thrustLR2[i])
    t2=sum(thrustLRstd2[i])
    thrust2.append(t1)
    thruststd2.append(t2)

#Cl and Cd for measurement series 2
cl2nd = list()
cd2nd = list()
for i,x,y,z in zip(thrust2,tas,dens,column(datadefl,8)):
    cl2nd.append(cl(weight(z),y,x))
    cd2nd.append(cd(i,y,x))
#print(cl2nd)
#print(cd2nd)



#Thrust coefficients
tc, tcs = list(),list()
for i,j,x,y in zip(thrust2,thruststd2,dens,redveq):
    tc.append(thrustcoeff(i,x,y,0.686))
    tcs.append(thrustcoeff(j,x,y,0.686))

####
###
### CG shift data
#getting data from CG shift measurements
CGshift = np.genfromtxt('CGShift.txt', dtype=float, comments='#',  delimiter="&", skip_header=1)

#Pressure
presscg = list()
for i in column(CGshift,0):
        presscg.append(pressure(i))
        
#mach number
machcg = list()
for i,x in zip(column(CGshift,1), presscg):
    machcg.append(Mach(x,i))

#temperature
tempercg = list()
for i,x in zip(column(CGshift,9), machcg):
    tempercg.append(temp(i,x))
 
#speed sound
ssoundcg = list()
for i in tempercg:
    ssoundcg.append(speedofsound(i))

#density
denscg = list()
for i,x in zip(presscg,tempercg):
    denscg.append(density(i,x))
    
#True Airspeed
tascg = list()
for i,x in zip(machcg,ssoundcg):
    tascg.append(TAS(i,x))

#CL
clcg = list()
for x,y,z in zip(tascg,denscg,column(CGshift,8)):
    clcg.append(cl(weight(z),y,x))

# Value of CN to be used to calculate Cmdelta
CN = sum(clcg)/2    #takes average of CL during cg shift experiment
changedefl = column(CGshift,3)[1]-column(CGshift,3)[0]
# Calculating Cmdelta    
cmd = Cmdelta(CN,deltaxcg,np.radians(changedefl)) # note radians for change in deflection

###
###
### reduced Elevator deflection
redeldef = list()
for x,y,z in zip(column(datadefl,3),tcs,tc):
    redeldef.append(redelevdef(x,y,z,cmd))
sortedredeldef = sorted(redeldef, key = float)    
# reduced stick force
redstickforce = list()
for x,y in zip(column(datadefl,5),column(datadefl,8)):
    redstickforce.append(redelvforce(x,y))
sortedredstickforce = sorted(redstickforce, key=float) 
###
###
### Plots and interpolation


# Set plot size
plt.rcParams.update({'font.size': 15})

#reduced equivalent velocity and reduced stickforce
xspace = np.linspace(65., 105., 1000)
coefs = poly.polyfit(sortedredveq, sortedredstickforce,3) # polynomial fit
ffit = poly.polyval(xspace, coefs)
plt.grid()
plt.ylabel('Reduced stick-force [N]')
plt.xlabel('Reduced equivalent airspeed [m/s]')
plt.scatter(sortedredveq, sortedredstickforce, label='measurement points')
plt.plot(xspace,ffit,"r--",label='polynomial fit')
plt.legend(loc=1)
plt.gca().invert_yaxis()
plt.show()

#reduced equivalent velocity and reduced elevator deflection
coefs1 = poly.polyfit(sortedredveq,sortedredeldef,3) # polynomial fit
ffit1 = poly.polyval(xspace, coefs1)
plt.grid()
plt.ylabel('Reduced elevetor deflection [deg]')
plt.xlabel('Reduced equivalent airspeed [m/s]')
plt.scatter(sortedredveq, sortedredeldef, label='measurement points')
plt.plot(xspace,ffit1,"r--",label='polynomial fit')
plt.legend(loc=1)
plt.gca().invert_yaxis()
plt.show()

#alpha and reduced elevator deflection 
xspace2 = np.linspace(2.8,8.2,100)
alpha = sorted(column(datadefl,2),key=float, reverse=True)
coefs2 = poly.polyfit(alpha,sortedredeldef,1) # polynomial fit - linear
ffit2 = poly.polyval(xspace2, coefs2)
plt.grid()
plt.ylabel('Reduced elevator deflection [deg]')
plt.xlabel('Angle of attack [deg]')
plt.scatter(alpha, sortedredeldef,label='measurement points')
plt.plot(xspace2,ffit2,"r--",label='polynomial fit')
plt.legend(loc=2)
plt.gca().invert_yaxis()
plt.show()

###
###
###
slope = poly.polyder(coefs2)  #slope of graph alpha vs reduced elevator deflection 
cmalpha = -cmd*slope # Cmalpha


cmdpercent = abs((cmd-cmdtheory)/cmdtheory)*100
cmalphapercent = abs((cmalpha-cmalphatheory)/cmalphatheory)*100
print('cmdelta = ',cmd)
print('cmdelta deviation vs theory ',cmdpercent,'%')
print('cmalpha = ',cmalpha)
print('cmalpha deviation vs theory ',cmalphapercent,'%')
