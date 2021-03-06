# -*- coding: utf-8 -*-
"""compareDISK.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1XFHByE-5kdhUnZYRfj2zL1k7NKV4U7Mz
"""

import matplotlib as mpl
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kn
from scipy.special import kvp
from scipy.interpolate import griddata
from matplotlib import cm
mpl.rcParams['legend.fontsize'] = 10

def extract_data_SE(filename):
  """
  Extracts the vertex data from a dump file created by
  surface evolver.
  """
  verticies_flag = False
  point_data = []
  parameters = []
  with open(filename,'r') as data:
    for line in data:
      try:  
        line = line.split()
        marker = line[0]
        
        # Find where verticies end
        if marker == "edges": 
          verticies_flag = False
          break

        # Extract data of verticies
        elif verticies_flag == True: 
          line = [float(x) for x in line[1:4]]
          #if line == [0,0,0]:  continue
          point_data.append(line)

        # Find where the verticies section starts
        elif marker == "vertices": 
          verticies_flag = True

        # Extract parameters
        elif marker == "PARAMETER":
          parameters.append(float(line[3]))
        elif marker == "GRAVITY_CONSTANT:":
          parameters.append(float(line[1]))

      except IndexError:
        continue

  point_data = np.array(point_data) # np arrays are better
  return point_data, parameters

def extract_data_FCD(filename):
  """
  Extract data from the dump of the FCD analysis 
  """
  point_data = []
  with open(filename,'r') as data:
    for line in data:
      try:  
        line = line.split()
        marker = line[0]
        
        line = [float(x) for x in line[1:4]]
        #if line == [0,0,0]:  continue
        point_data.append(line)

      except IndexError:
        continue

  point_data = np.array(point_data) # np arrays are better
  return point_data, parameters

def z_Harry(r,phi,alpha,l_c,R_in):
  lratio = R_in/l_c
  Q = lratio/8 * kvp(1,lratio) / kn(1,lratio)
  A = 1
  B = 3*Q - 1/24
  D = Q - 1/8

  z = R_in*(alpha)* ((A+B*alpha**2)*(kn(1,r/l_c)/kn(1,lratio)) * np.cos(phi)+
                         ((D*alpha**2)*(kn(3,r/l_c)/kn(3,lratio)) * np.cos(3*phi)))
  
  return np.nan_to_num(z,posinf=0.0,neginf=0.0) # remove inf errors
  
def z_Dom(r,phi,alpha,l_c,R_in):
  lratio = R_in/l_c
  z = (alpha)* R_in *(kn(1,r/l_c)/kn(1,lratio))* np.cos(phi) # Solution
  return np.nan_to_num(z,posinf=0.0,neginf=0.0) # remove inf errors

def z_Experimental(r,phi,alpha,l_c,R_in):
  lratio = R_in*np.cos(alpha)/l_c
  z = np.sim(alpha)* R_in *(kn(1,r/l_c)/kn(1,lratio))* np.cos(phi) # Solution
  return np.nan_to_num(z,posinf=0.0,neginf=0.0) # remove inf errors

def absDif(z1,z2,z_max):
  return abs((z1-z2)/z_max)*100

from numpy.lib.type_check import mintypecode
from numpy.core.arrayprint import DatetimeFormat
point_data, parameters = extract_data_SE('/content/drive/MyDrive/Part III Project/Surface Evolver/alpha/0.065.txt')
x, y, z_SE = point_data.T

R_in, R_out, alpha, gamma, rho, g = parameters

l_c = (gamma/(rho*g))**0.5 # Capillary length sqrt(gamma/rho g)

r = (x**2 + y**2)**0.5
phi = np.arctan2(y,x)

z_h = z_Harry(r,phi,alpha,l_c,R_in)
z_d = z_Dom(r,phi,alpha,l_c,R_in)
z_max = R_in * np.sin(alpha)

print('Expected max:', z_max)
print('z max:',z_SE.max(), ' dif:', z_max-z_SE.max())
print('z_dom max:', z_d.max(), ' dif:', z_max-z_d.max());
print('z_harry max:', z_h.max(), ' dif:', z_max-z_h.max());
print(l_c)

fig, ax = plt.subplots(figsize=(10,10))
q= ax.tripcolor(x, y, absDif(z_SE,z_d,z_max))
print(max(absDif(z_SE,z_d,z_max)))
cbar=fig.colorbar(q)
ax.axis('scaled')
#ax.set(xlim=(R_in-0.01, R_in+0.01), ylim=(-0.01, 0.01));
ax.set(xlim=(-R_out, R_out), ylim=(-R_out, R_out));

fig, ax = plt.subplots(figsize=(10,10))
q= ax.tripcolor(x, y, absDif(z_SE,z_h,z_max))
print(max(absDif(z_SE,z_h,z_max)))
cbar=fig.colorbar(q)
ax.axis('scaled')
#ax.set(xlim=(R_in-0.01, R_in+0.01), ylim=(-0.01, 0.01));
ax.set(xlim=(-R_out, R_out), ylim=(-R_out, R_out));

"""
kept alpha const

"""
# N's tested
N_tested = [3,4,5,8,25,12,16,20,10,40]

# Max error taken from above graph
error_tested_dom = [15.258678528493224,10.919861120893652,7.742950709890792,2.6835077612798597,
                    1.5218051000763697,0.6344387750319269,0.7840802630950636,1.0866916246477476,
                    1.3080861519563247,2.024273347521723]

# Max error taken from above graph
error_tested_harry = [15.236153350978087,10.901844206728796,7.728749862815447,2.676999958303302,
                      1.148854374319932,0.6323683513970342,0.471205865898834,0.7432962318694727,
                      1.304372155326567,1.5732245782659735]
error_tested_longd = [14.64389799609555,0,1.0625692395381938,0.6735310007991866,0,0.39632152467025333,
                      0.5322828109692251,0.6795764772030947,0.3521525530468352,0]
error_tested_longh = [14.732654620995705,0,1.1452530269384882,0.7802948113165578,0,0.10824171356102062,
                      0.3439020374145709,0.28060314002016923,0.31751919422498376,0]

m, A = np.polyfit(np.log(N_tested), np.log(error_tested_harry), 1)
print('Fitted exponent:', m,A)

fig, ax = plt.subplots()

ax.loglog(N_tested, error_tested_harry, 'r+',ms=10)
ax.loglog(N_tested, error_tested_dom, 'k+',ms=10)
ax.loglog(N_tested, error_tested_longd, 'g+',ms=10)
ax.loglog(N_tested, error_tested_longh, 'b+',ms=10)

ax.set_ylabel('Maximum error')
ax.set_xlabel('R_in')
ax.legend(['harry','dom','longd','longh'])

plt.show();

"""
kept R_max constant and varied R/alpha to plot max errors
clearly larger r are better but htis may be due to the smaller angles

"""
# N's tested
N_tested = [3,4,5,8,12,16,20,25,30,50,100]

# Max error taken from above graph
error_tested_dom = [58.64379069310925,28.60897131210795,18.226671471857028,8.194621455149068,4.533951260988207,3.1954033327969547,2.6123953831696265,2.4182980115678165,2.1722144954826508,1.9834411121191071,4.275401125402529]

# Max error taken from above graph
error_tested_harry = [28.152468639075035,8.445727212171503,10.118564826108447,1.3295626497814628,1.15950099327983,1.2966104956968176,1.3410021605501807,1.2136123148148412,1.2221925162793348,1.5647598640462335,4.245353489770937]

m, A = np.polyfit(np.log(N_tested), np.log(error_tested_harry), 1)
print('Fitted exponent:', m,A)

fig, ax = plt.subplots()

ax.loglog(N_tested, error_tested_harry, 'r+',ms=10)
ax.loglog(N_tested, error_tested_dom, 'k+',ms=10)

ax.set_ylabel('Maximum error')
ax.set_xlabel('R_in')
ax.legend(['harry','dom'])

plt.show();

print(l_c)
error_h = 10e10
error_d = 10e10
l_c_min_h = 0
l_c_min_d = 0
for l_cc in  np.linspace(l_c-0.8*l_c,l_c+0.8*l_c,200):
  z_h_t = z_Harry(r,phi,alpha,l_cc,R_in)
  z_d_t = z_Dom(r,phi,alpha,l_cc,R_in)
  
  z_dif = absDif(z_SE,z_h_t,z_max)
  if sum(z_dif**2) < error_h:
    error_h = sum(z_dif**2)
    l_c_min_h = l_cc
  z_dif = absDif(z_SE,z_d_t,z_max)
  if sum(z_dif**2) < error_d:
    error_d = sum(z_dif**2)
    l_c_min_d = l_cc

print('Dom',l_c_min_d, abs(l_c-l_c_min_d)/l_c*100)
print('Harry',l_c_min_h, abs(l_c-l_c_min_h)/l_c*100)

rs = np.linspace(R_in, R_in+0.01, 100)
ps = 0#np.linspace(0, 2*np.pi, 100) # Angles
R, P = np.meshgrid(rs, ps)

grid_se = griddata(np.array([x,y]).T, z_SE, (R*np.cos(P), R*np.sin(P)), method='cubic')

grid_harry_min = griddata(np.array([x,y]).T, z_Harry(r,phi,alpha,l_c_min_h,R_in), (R*np.cos(P), R*np.sin(P)), method='cubic')
grid_harry = griddata(np.array([x,y]).T, z_h, (R*np.cos(P), R*np.sin(P)), method='cubic')

grid_dom_min = griddata(np.array([x,y]).T, z_Dom(r,phi,alpha,l_c_min_d,R_in), (R*np.cos(P), R*np.sin(P)), method='cubic')
grid_dom = griddata(np.array([x,y]).T, z_d, (R*np.cos(P), R*np.sin(P)), method='cubic')

grid_exp = griddata(np.array([x,y]).T, z_Dom(r,phi,alpha,l_c,R_in), (R*np.cos(P), R*np.sin(P)), method='cubic')

fig, ax = plt.subplots(figsize=(20,10))

ax.plot(1000*(rs-R_in),grid_se[0])
ax.plot(1000*(rs-R_in),grid_dom[0],'.')
ax.plot(1000*(rs-R_in),grid_dom_min[0],'-.')
ax.plot(1000*(rs-R_in),grid_harry[0],'.')
ax.plot(1000*(rs-R_in),grid_harry_min[0],'-.')
ax.plot(1000*(rs-R_in),grid_exp[0],'-.')

ax.set_xlabel('r-R_in / mm')
ax.set_ylabel('h / mm')
ax.legend(['SE','Dom','Dom_min','Harry','Harry_min','exp'])
plt.yscale("log")
plt.show()

# square sum
# N's tested
N_tested = [0.1,0.05,0.08,0.15,0.065,0.12]

# Max error taken from above graph
error_tested_dom = [14.07035175879397,3.618090452261312,9.246231155778903,30.76923076923077,6.0301507537688375,19.698492462311563]

# Max error taken from above graph
error_tested_harry = [10.050251256281397,2.8140703517587875,6.834170854271362,22.564102564102555,4.4221105527638205,14.87437185929648]

m, A = np.polyfit(np.log(N_tested), np.log(error_tested_harry), 1)
print('Fitted exponent:', m,A)

fig, ax = plt.subplots()

ax.loglog(N_tested, error_tested_harry, 'r+',ms=10)
ax.loglog(N_tested, error_tested_dom, 'k+',ms=10)

ax.set_ylabel('Maximum error')
ax.set_xlabel('alpha')
ax.legend(['harry','dom'])

plt.show();

# N's tested
N_tested = [0.1,0.2,0.05,0.08,0.15,0.065,0.12]

# Max error taken from above graph
error_tested_dom = [0.10413218763468281,0.6970363176096418,0.026382672706357333,0.06827167214908429,0.2824921003726683,0.044081212239643244,0.15867563318246566]

# Max error taken from above graph
error_tested_harry = [0.06954247686555401,0.44851839079366784,0.020291344685508104,0.04565517560903049,0.18451204618337386,0.030997856858919096,0.1058623701854268]

m, A = np.polyfit(np.log(N_tested), np.log(error_tested_harry), 1)
print('Fitted exponent:', m,A)

fig, ax = plt.subplots()

ax.loglog(N_tested, error_tested_harry, 'r+',ms=10)
ax.loglog(N_tested, error_tested_dom, 'k+',ms=10)

ax.set_ylabel('Maximum error')
ax.set_xlabel('alpha')
ax.legend(['harry','dom'])

plt.show();
