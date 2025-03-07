import csv
from rich import print
import math
import numpy as np
import itertools
import matplotlib.pyplot as plt

kpctocm = 3.086e21 ## convert kpc to cm
pc = kpctocm*1e-3 ## convert pc to cm 
Msun = 1.989e33   ## in g
## location of Sun -
rsun= 8.122 ## in kpc
xsun = rsun
ysun = 0.
zsun = 0.0055  ## in kpc from Quillen et al. 2020
phisun = np.arctan2(ysun, xsun) ## should be in radian 
xsun *= kpctocm
ysun *= kpctocm
zsun *= kpctocm
rsun *= kpctocm

G = 6.67e-8 # in cgs
c = 3e10  ## in cm
Vlsr = 255.2*1e5 ## in cm/s
alpha1_0 = 4e-30
lgalpha1 = math.log10( alpha1_0)

def get_al(item):
    x, y, z = findxyz(item)
    measured_alos = float(item["ALOS_PS"])*(1./(10.*3.14e7))  ## this is a single number in cm/s^2
    alos_Galactic = alos_gal(x, y, z) ## this is the model Galactic acceleration
    alos_sun = alos_gal(xsun,ysun,zsun)

    dx, dy, dz = x-xsun, y-ysun, z-zsun
    dr = np.sqrt(dx*dx+dy*dy+dz*dz)
    ax,ay,az = tuple(a - b for a, b in zip(alos_Galactic, alos_sun))

    al = (ax*dx + ay*dy + az*dz)/(np.maximum(dr,1.e-10)) ## this is the Galactic LOS acceleration (with the Sun's acceleration subtracted out) at vector X1.

    return al

def alos_gal(x,y,z):
    ## assume VERY SIMPLE Quillen model for alos_gal 
    r = np.sqrt(x*x + y*y)
    alpha1 = 1e1**lgalpha1
    alpha2 = 0.

    az = -alpha1*z - alpha2*(z*z)*np.sign(z)
    ar = Vlsr*Vlsr/r

    ax = -ar*x/r
    ay = -ar*y/r
    #print ("x,y,z",x,y,z)
    
    return ax, ay, az 

def radians(degrees):
    return degrees * (math.pi / 180)

def findxyz(item):
    l_rad = radians(float(item["GL"]))
    b_rad = radians(float(item["GB"]))

    # Convert parallax from mas to arcseconds
    p_arcsec = float(item["PX"]) / 1000.0

    # Calculate the distance in parsecs
    distance = 1 / p_arcsec  # distance in parsecs

    x = distance * math.cos(b_rad) * -math.cos(l_rad)
    y = distance * math.cos(b_rad) * math.sin(l_rad)
    z = distance * math.sin(b_rad)
    x = x*pc + xsun # now in cm 
    y = y*pc + ysun
    z = z*pc + zsun 

    return x,y,z

# This csv file has the spin-period measurements. ALOS_PS has the Schlovskii and magnetic braking terms subtracted out.  Units are mm/s/yr = cm/s/decade.  From: https://github.com/thomasdonlon/Empirical_Model_MSP_Spindown_Accels/blob/main/data.csv
with open("data.csv", mode="r") as file:
    # Create a DictReader object
    csv_reader = csv.DictReader(file)

    # Iterate over each row in the CSV
    items = []
    for row in csv_reader:
        #item = {k: float(v) for k, v in row.items() if k in ["GL", "GB", "PX", "ALOS_PS"]}
        ## ignore the empty rows for alos_ps -
        #item = {k: float(v) for k, v in row.items() if k in ["GL", "GB", "PX", "ALOS_PS"] and v != ''}
        ## set alos_ps to 0. when it is not present (otherwise problem w above):
        #item = {k: float(v) if v != '' else 0. for k, v in row.items() if k in ["GL", "GB", "PX", "ALOS_PS"]}
        item = {k: v if v != '' else 0. for k, v in row.items() if k in ["GL", "GB", "PX", "ALOS_PS", "NAME"]}

        items.append(item)
        
