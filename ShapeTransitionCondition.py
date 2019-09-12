#PROGRAM TO CALCULATE NECESSARY CONDITION ON K24 TO ALLOW FOR A SHAPE TRANSITION.
#WRITTEN TO ACCOMPANY "MOLECULAR HETEROGENEITY INDUCES RECONFIGURABLE NEMATIC LIQUID CRYSTAL DROPS"
#WRITTEN BY SOPHIE ETTINGER AND WEI-SHAO WEI SEPT 2019.

from  __future__ import print_function
import numpy as np
from scipy import optimize

#FREE ENERGY EQUATIONS USING ONE ELASTIC CONSTANT
def F1(a):
	K24 = a*K
	sigma = (Wa*r/K)+(K24/K)-1.
	#F_s IS FREE ENERGY OF A SPHERE ASSUMING RADIAL CONFIGURATION
	F_s = 8.*np.pi*R*(K-0.5*K24)+gamma*4.*np.pi*R**2
	#F_f IS FREE ENERGY OF A FILAMENT ASSUMING ESCAPED RADIAL CONFIGURATION
	F_f = np.pi*K*L*(3.-(K24/K)-(1./sigma))+gamma*2.*np.pi*r*L
	#SET TWO EQUATIONS EQUAL TO FIND WHERE TRANSITION ENERGETICALLY OCCURS
	return F_f-F_s

#FREE ENERGY EQUATIONS USING THREE ELASTIC CONSTANTS (K11, K33, & K24)
def F3(a):
	K24 = a*K
	sigma = (Wa*r/K)+(K24/K)-1.
	#F_s IS FREE ENERGY OF A SPHERE ASSUMING RADIAL CONFIGURATION
	F_s = 8.*np.pi*R*(K-0.5*K24)+gamma*4.*np.pi*R**2
	#F_f IS FREE ENERGY OF A FILAMENT ASSUMING ESCAPED RADIAL CONFIGURATION
	if k > 1:
		k0 = np.sqrt(k-1.)
		F_f = (gamma*2.*np.pi*r*L)+(np.pi*K*L*(2.+((k/k0)*np.arctan(k0))-((k/k0)*np.arctan(k0/sigma))-(K24/K)))
	if k < 1:
		k0 = np.sqrt(1.-k)
		F_f = (gamma*2.*np.pi*r*L)+(np.pi*K*L*(2.+((k/k0)*np.arctanh(k0))-((k/k0)*np.arctanh(k0/sigma))-(K24/K)))
	#SET TWO EQUATIONS EQUAL TO FIND WHERE TRANSITION ENERGETICALLY OCCURS
	return F_f-F_s

#ASSUME REASONABLE VALUES FOR R AND r (IN METERS)
#R AND r HERE ARE TAKEN FROM REAL EXPERIMENTAL DATA
R = 13.e-6
r = 0.65e-6
#CALCULATE FILAMENT LENGTH ASSUMING CONSTANT VOLUME
L = (4./3.)*(R**3)/(r**2)
#INITIALIZE CONDITIONS
cond1 = False
cond2 = False


#PROGRAM BEGINS HERE
print (' ')
print ('This program will calculate the shape transition condition for any given set of K (elastic constant), Wa (anchoring energy coefficient), and gamma (interfacial tension).')
print (' ')
print ('For detailed free energy equations, consult Methods section of the paper, Eqn. (1)-(7).')
print (' ')
print ('Note, consult Methods section of the paper for suggested values of parameters. e.g., K=1e-10; Wa=1e-4; gamma=0.236e-3')
i = int(input('Enter 1 to use the one-constant approximation or 3 to use all elastic constants: '))
print ('Note: Please enter values as a decimal or in the form of 1e-10, without spaces.')

#USE K OR K11 DEPENDING ON WHICH APPROXIMATION WE USE
if i == 1:
	Kstring = 'K'
if i ==3:
	Kstring = 'K11'

#USER-DEFINED VALUES FOR ELASTIC CONSTANTS, Wa, AND GAMMA. LOOP UNTIL BOTH CONDITIONS ARE MET.
while cond1==False or cond2==False:
	K = float(input("Enter value for %s, in N: " % Kstring))
	if i ==3:
		k = float(input("Enter value for ratio of K33/K11: "))
	Wa = float(input("Enter value for Wa, in J/m^2: "))
	gamma = float(input("Enter value for gamma, in N/m: "))
	if K/Wa <= 2.6e-6:
		cond1 = True
	else:
		print ("%s/Wa too large; choose a different combination." % Kstring)
	if (gamma*R)/(2.*K) >= 10.:
		cond2 = True
	else:
		print ("gamma/%s too small; choose a different combination." % Kstring)

print (' ')

#SOLVE FOR SHAPE TRANSITION CONDITION USING ONE CONSTANT APPROXIMATION
#CHOOSE 5 AT INITIAL GUESS.
if i == 1:
	a_value = optimize.fsolve(F1,5)
	print ('If K24 > %.2f * K, shape transition is possible.' % a_value)

#SOLVE FOR SHAPE TRANSITION CONDITION USING THREE CONSTANTS
#CHOOSE 5 AS INITIAL GUESS.
if i == 3:
	a_value = optimize.fsolve(F3,5)
	print ('If K24 > %.2f * K11, shape transition is possible.' % a_value)

