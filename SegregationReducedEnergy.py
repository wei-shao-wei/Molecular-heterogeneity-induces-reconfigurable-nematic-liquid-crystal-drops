#PROGRAM TO CALCULATE ELASTIC FREE ENERGY DIFFERENCE (REDUCTION) DUE TO OLIGOMER SEGREGATION.
#WRITTEN TO ACCOMPANY "MOLECULAR HETEROGENEITY INDUCES RECONFIGURABLE NEMATIC LIQUID CRYSTAL DROPS"
#WRITTEN BY SOPHIE ETTINGER AND WEI-SHAO WEI SEPT 2019.

from  __future__ import print_function
import numpy as np
from scipy import optimize

#THE MOLECULAR WEIGHT OF THE MACROMER IS 6641 g/mol.
#THE MOLECULAR WEIGHT OF THE MONOMER IS 673 g/mol.
#f_s=V_shell/V_total IS THE VOLUME FRACTION OF THE SHELL, 
#f_c=V_core/V_total IS THE VOLUME FRACTION OF THE CORE; SUCH THAT f_s+f_c=1.
def func(n_val):
	#A IS THE WEIGHT OF THE MACROMER IN THE SHELL
	A = (f_s*6641.*(l_s-1.))/(6641.*(l_s-1.)+673.*(9.-l_s))
	#B IS THE WEIGHT OF THE MACROMER IN THE CORE
	B = (f_c*6641.*n_val)/((6641.*n_val)+(673.*(1.-n_val)))
	#C IS THE WEIGHT OF THE MACROMER IN THE ENTIRE STRUCTURE
	C = w_macro/(1+w_macro)
	#THE TOTAL AMOUNT OF WEIGHT OF MONOMER AND MACROMER IS CONSERVED DURING SEGREGATION, i.e., C=A+B
	return A+B-C

#PROGRAM BEGINS HERE
#THE AVERAGE CHAIN LENGTH FOR THE MONOMER IS <l> = 1 AND THE AVERAGE CHAIN LENGTH FOR THE MACROMER IS <l> = 9.
#FOR THESE CALCULATIONS, WE ASSUME K IS PROPORTIONAL TO THE AVERAGE CHAIN LENGTH.

print (' ')
print ('This code will show how the macromers and monomers segregate into core and shell regions with uneven chain length distributions in order to lower the system\'s overall elastic energy. Here we assume K is proportional to <l>.')
print (' ')
print ('Please consult Methods section of the paper for suggested values of the parameters. e.g., macromer:monomer weight ratio=0.7; mean chain length of shell=1.7; volume ratio V_shell/V_total=0.5.')
print (' ')
print ('Note: Please enter all values as decimals.')

#USER-DEFINED VALUES FOR WEIGHT RATIO OF MACROMER IN THE RANGE [0.7,1.0], BASED ON THE AVAILABLE EXPERIMENTS EXHIBITING FILAMENTS.
w_macro = float(input('Enter the weight ratio of macromer:monomer (choose a number between 0.7 and 1.0): '))

#CALCULATE AVERAGE MEAN CHAIN LENGTH w/ GIVEN WEIGHT RATIO OF MACROMER:MONOMER
#THE MOLECULAR WEIGHT OF THE MACROMER IS 6641 g/mol.
#THE MOLECULAR WEIGHT OF THE MONOMER IS 673 g/mol.
l_avg=(1.*(1./673.)+9.*(w_macro/6641.))/(1./673.+w_macro/6641.)

print ('We start with a monomer:macromer weight ratio of 1:%.1f. This gives an overall mean chain length of <l> = %.2f.' % (w_macro,l_avg))
print (' ')
print ('Note, we employ the chain length <l> = 1 for the monomer and <l> = 9 for the macromer.')

#USER-DEFINED VALUES FOR THE MEAN CHAIN LENGTH IN THE LONG-CHAIN RICH SHELL, l_s, AND THE VOLUME FRACTION OF THE SHELL, f_s. 
l_s = float(input('Enter the mean chain length, <l>, of the long-chain-rich shell: '))
f_s = float(input('Enter ratio of shell volume to the total drop volume (V_shell/V_total): '))

#FIND f_c, THE VOLUME FRACTION OF THE SHORT-CHAIN RICH CORE.
f_c = 1-f_s

#n IS THE NUMBER FRACTION OF MACROMERS IN THE SHORT-CHAIN RICH CORE REGION.
#FIRST WE FIND THE ROOTS OF THE DEFINED FUNCTION.
print (' ')
n = optimize.fsolve(func,1)
l_c = (1.-n)+9.*n
print ('The mean chain length of the short-chain-rich core is <l> = %.2f' % l_c)
print (' ')

#FOR A SPHERE: THIS IS THE ELASTIC FREE ENERGY RATIO OF SEGREGATED:HOMOGENEOUS. HERE WE USE THE FREE ENERGY FORM 8*Pi*K*R.
#WE IGNORE THE K24 TERM.
Fs = ((8.*n+1.)*np.cbrt(f_c)+l_s*(1.-np.cbrt(f_c)))/l_avg
print ('For a spherical drop with a radial configuration, the segregated case will have an elastic energy that is %.1f%% of the homogenous case.' % (Fs*100.))
print (' ')

#FOR A CYLINDRICAL FILAMENT WITH PLANAR RADIAL CONFIGURATION: THIS IS THE ELASTIC FREE ENERGY RATIO OF SEGREGATED:HOMOGENEOUS. HERE WE USE THE FREE ENERGY FORM Pi*K*ln(r/rho), WHERE RHO IS THE SIZE OF THE DEFECT CORE.
#USUALLY THE CORE SIZE IS AROUND 10nm AND OUR FILAMENTS HAVE AN AVERAGE RADIUS OF 1um, SO WE FIND A RATIO r/rho OF ~100. 
#WE IGNORE THE K24 TERM.
#NOTE THAT WE USE THE PLANAR RADIAL CONFIGURATION TO APPROXIMATE THE ESCAPE RADIAL CONFIGURATION.
f_c0 = np.sqrt(f_c)
Ff = ((8.*n+1.)*np.log(100.*f_c0)+l_s*np.log(1./f_c0))/(l_avg*np.log(100.))
print ('For a cylindrical filament with planar radial (not escape radial) configuration and a defect core size of ~10nm, the segregated case will have an elastic energy that is %.1f%% of the homogenous case.' % (Ff*100.))


