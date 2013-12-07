'''
	lens.py

Uses the fortran backend lc_sub to calculate the magnification
for a binary lens which is prepared for python using f2py

e.g. f2py -c lc_sub.f -m lc_sub

'''

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import lc_sub
import time

# define a few constants
G = 6.67e-11
C = 299792458.0
PC = 3.086e16
AU = 1.5e11
MSUN = 1.98892e30

# here you specify the values of
# v velocity in m/s, the two masses m2 and m1 in solar masses, the distance d_l in pc to the lens
# the two inclination angles, omega and i_p (i_p=0.5pi is edge on), 
# and the value of beta, clock specifies anticlockwise or clockwise orbit

def do_light_curve ( m1, m2, alpha, beta, d_l, v, i_p, omega, dt = 1.0, tmax=48000.0, phi_0 = 0.5):
	'''calculate a lensing light curve'''

	mtot = m1 + m2
	a1 = alpha * m2 / mtot
	a2 = alpha * m1 / mtot

	# theta_e in arcseconds
	theta_e = 2.0 * np.sqrt((G*mtot*MSUN)/(d_l * PC* C**2)) * (180.0/np.pi) * 3600.0

	# Einstein radius and semi major axis in metres
	R_E = theta_e * d_l
	a = alpha * R_E    

	# orbital period, days             
	P = 365.25 * np.sqrt (a**3 / mtot)
      
	# orbital period, hours      
	p_hour = P * 24.0

	# proper motion and initial position
	mu = ( ( v / AU) * (3600.0) )  /  (d_l * theta_e)
	x_i = - ( mu * 24000.0 ) + 0.000000001

	phi_0 *= np.pi

	time = np.arange ( 0.0, tmax, dt )
	w = 2.0 * np.pi / p_hour

	nrnew = 3

	mags = []
	pos1 = []
	pos2 = []

	for t in time:

		phi = phi_0 - (w * t) 
		xcm = x_i + (mu * t)
		ycm = beta
          
		# now work out the position of mass 1 (x1, y1) and mass 2 (x2, y2)
		# these are a function of time (phase: phi) and orientation (i_p and omega)       
		x1 = xcm - a1 * (  ( np.cos(omega)*np.cos(phi) ) - \
				( np.sin(omega)*np.sin(phi)*np.cos(i_p) )  )
		y1 = ycm - a1 * (  ( np.sin(omega)*np.cos(phi) ) + \
				( np.cos(omega) * np.sin(phi) * np.cos(i_p) )  )
		x2 = xcm + a2 * (  ( np.cos(omega)*np.cos(phi) ) - \
				( np.sin(omega)*np.sin(phi)*np.cos(i_p) )  )
		y2 = ycm + a2 * (  ( np.sin(omega)*np.cos(phi) ) + \
                ( np.cos(omega)*np.sin(phi)*np.cos(i_p) )  )

		pos1 .append ( [x1, y1] )
		pos2 .append ( [x2, y2] )

		u_cm = np.sqrt ( xcm**2 + ycm**2)

		# magnification of hypothetical single lens          
		A_s = mag_single (u_cm) 

		# source coordinates is origin          
		x = 0.0
		y = 0.0 
                    
		A_b = lc_sub.mag( x, y, x1, y1, x2, y2, m1/mtot, m2/mtot, nrnew)

		mags.append (A_b)

	mags = np.array (mags)

	return time, mags

		
	


def plot_light_curve ( m1, m2, alpha, beta, d_l, mu, i_p, omega, dt = 1.0, tmax=48000.0, phi_0 = 0.5):
	''' plot light curve up'''
	
	import matplotlib.pyplot as plt

	time, mag = do_light_curve ( m1, m2, alpha, beta, d_l, v, i_p, omega, dt = dt, tmax= tmax, phi_0 = phi_0)

	fig = plt.figure()

	ax = fig.add_subplot(111)
	ax.plot(time, mag)
	ax.set_xlabel ("Time (hours)")
	ax.set_ylabel ("Magnification")

	plt.show()



def mag_single ( u ):
	'''magnification of single lens'''
	A = ( u**2 + 2.0) / ( u * np.sqrt( u**2 + 4.0) )
	return A

n_lightcurves = 10


for i in range( n_lightcurves):   
   
	m2=1.3 * np.random.random()
	m1=0.3 * np.random.random()
	alpha=0.35 * np.random.random()
	d_l=50.0 * np.random.random()
	v=10000.0 * np.random.random()
		
	# These next two values determination the orientation of the system
	# i_p being inclination and omega being ascending node longitude
	# consult figure 
	i_p = np.pi * np.random.random() 
	omega = np.pi * np.random.random()
		  
	# beta is the angular distance of closest approach      
	beta=0.5

	print "Lightcurve %i" % i
	time1 = time.time()
	do_light_curve ( m1, m2, alpha, beta, d_l, v, i_p, omega)
	time2 = time.time()
	print "Time %.2f seconds" % (time2 - time1)





