c-----------------------------------------------------------------------------------
c
c	    James Matthews, Rosanne Di Stefano, Shude Mao
c	    March 26th 2013, Southampton University (based on previous codes)
c
c			------SPIDER.f------
c
c	This code is intended to be an all encompassing fortran code for 
c	elliptical and inclined orbits. We compute light curves for a large 
c	number of different types of orbits, and the code identifies 
c	maxima in the deviation from the standard point lens light curve.
c
c	Step 1: Specify general parameters, e.g masses, lens distance, PM, beta
c	Step 2: Iterate over randonly generated phase and orbital elements
c	Step 3: Compute light curve along path, and identify peaks
c	Step 4: Write peaks and light curves to file
c
c	NOTE: THIS DOES NOT WORK FOR ECCENTRICITY YET!!
c-----------------------------------------------------------------------------------

c	Define and assign variables- some of these variables are a bit messy.

	Implicit none

c	parameters and constants:

	Real*8 G, msun, c, pi, pc, au, VERYBIG, SMALL
	Parameter(pi=3.14159265, G=6.67d-11, c=299792458.0d0, msun=1.98892d30, 
     *	          pc=3.086d16, au=1.5d11, VERYBIG=1.0d10, SMALL=1.0d-10) 

c	orbital variables, elliptical angles etc:

	Real*8 omega, i_p, phi, phi_0, true_anom, mean_anom, ecc_anom, e, ecc
	Integer i, it, iter, iteration, i_phase, i_inc, i_omega, idum, nr, nrnew, nr0, number_of_iterations

c	variables for RNG:

	Real xx, ran1	

c	system characteristics:

	Real*8 alpha, beta, a, a1, a2, m1, m2, mtot, d_l, d_s, theta, theta_e, r_e, einstein, P, p_day, ucm
	Real*8 x, y, z, x1, x2, y1, y2, xcm, ycm, u_pl, u_cm, A_b, A_s, mag, mag_s, aa, xdat, ydat
	Real*8 n_m, q, dev, peakdev, peakmag, u_peak, x_peak, y_peak, phi_peak, theta_e_star

c	arrays of positions and times:

	Integer itmax, tmax, sample, iter_max, max_iphase, max_iinc
	Parameter(tmax=17532)
	Real*8 times(tmax), xpos(tmax), years(tmax), probs(tmax), magpt(tmax), trues(tmax), mean(tmax), eccs(tmax)
	real*8 xstar(tmax), ystar(tmax), ypos(tmax), xcms(tmax), ycms(tmax), radius(tmax), days(tmax)
	character(3) pathlabel
	character(15) lc_filename, peak_filename
	CHARACTER(len=32) :: arg	
	Real*8 mu, tau_e, p_orb, w_peri, b1, b2, b, mu_teta_hr, x_init, r_e_au
	Real*8 ep, FF, Fd, EE, amp, time, t_max, u 
	Real*8 mag_dif,mag_dif_max,t_dif,phi_dif_mag_max, ref_now, ref_prev
	Real*8 t_late, t_early, u2, u1 , mag_1, mag_2
        Real*8 magpt_max, t_pt_max, deviation, dev_ref, dev_abs, tdev
	Real*8 mag_pt(tmax),mag_pt1,mag_pt2,mag_bin, dev1, dev2, dev_bin 
	Integer iphase,ialpha, iarr, imass, val0001, peak0001, caust, caust1 

c	these are the arrays of integers

	Integer t_pre(7), t_post(7), up(7), down(7), counter_peak, counter_vall,magup(7), magdown(7) 
	Real*8 alpha_max,alpha_min,p_max,p_min,mag_max,phi_mag_max, refs(7), trefs(7) 

c------------------------------------------------------------------------------------

c	print statement for header on screen.
c	you can edit this to write to an output file, for example (calls subroutine infowrite())

	call infowrite()

c	test routine for getting arguments and writing help
	DO i = 1, iargc()
          CALL getarg(i, arg)
c          WRITE (*,*) arg
	  if (arg.eq.'h'.or.arg.eq.'help') then
	    call help()
	  endif
        END DO

c	number of iterations you want to do
	number_of_iterations=20000


c------------------------------------------------------------------------------------

c       initialize the random number generator

	write(*,*) 'Initialising random number generator...'
        idum = -1
        do i=1, 75181
           xx = ran1(idum)
        enddo

c------------------------------------------------------------------------------------

c	open files for writing

	write(*,*) 'Opening files and assigning global parameters (masses, distances etc.)...'
	lc_filename='lc_inc_j'
	peak_filename='peaks_inc_j'
	open(31,file=lc_filename)
	Open(41,file=peak_filename)
	Open(10,file='test')

c------------------------------------------------------------------------------------

c	set up arrays an calculations
c	sample is the number of points we sample (i.e. every 1 points in this case)
c	we run 20000,50000 light curves (say) 
c	everytime there is a peak in the deviation from single lens behaviour we note the characteristics of that peak
c	for each light curve we record the number of peaks in given places, and the times associated with the events
c	we always note the number of images, nr or nrnew, so we can identify caustic crossings or 2-image dropouts

	sample=1	
	magpt_max=0.0
	beta=3.0d0	
	m2=0.075d0
	q=0.001/m2
	m1=q*m2
	mtot=m1+m2
	d_l=5.82

c	theta_e in arcseconds. neglect d_l/d_s

	theta_e=2.0d0*dsqrt((G*mtot*msun)/(d_l*pc*c**2))*(180.0d0/pi)*3600.0
	theta_e_star=2.0d0*dsqrt((G*m2*msun)/(d_l*pc*c**2))*(180.0d0/pi)*3600.0
	r_e_au=theta_e*d_l	
	r_e=r_e_au*au	

c	Let us step every hour. we'll do it for 2 years. 
c	Thats 17532 steps. mu is arcsec/year, x_init is 
c	starting position for COM ( 1 year from close approach) 
	          
	mu=1.5d0
	mu_teta_hr=(mu/(365.25*24.0))/theta_e
	x_init=-(mu/theta_e)
c	create reference arrays and explicitly set other to 0. ref is an array of reference points for the magnitudes. trefs is the same for the times.
	Do 501 iarr=1,7
	  refs(iarr)=0.01*(2.5119**real(iarr))
	  trefs(iarr)=2**real(iarr)
	  up(iarr)=0
	  down(iarr)=0
	  t_pre(iarr)=0
	  t_post(iarr)=0
	  magup(iarr)=0
	  magdown(iarr)=0
501	enddo
	Do 502 iarr=1,tmax
	  days(iarr)=real(iarr)/24.0
	  xcms(iarr)=x_init+mu_teta_hr*real(iarr)
	  u=sqrt(xcms(iarr)**2+beta**2)
	  mag_pt(iarr)=amp(u)
	  times(iarr)=real(iarr)
502	enddo

c------------------------------------------------------------------------------------

c	First we must define the orbital elements. 
c	with i_p and omega we can control whether the orbit is counter clockwise or not.

	imass=1
	write(*,*) 'Computing ', number_of_iterations,' Light Curves...'

	Do 100 i_inc=1,number_of_iterations

c	    we choose alpha randomly in log space, between 0.05, 20
	    alpha=10.0d0**(-1.301d0+(2.602d0*ran1(idum)))
	    a=r_e_au*alpha
	    p_day=365.25*dsqrt(a**3/mtot)
	    P=p_day*24.0d0

c	    omega is the longitude of ascending node (the angle that defines axis for inclination)
	    omega=2.0*pi*ran1(idum)

c	    i_p is the inclination
	    i_p=2.0*pi*ran1(idum)

c	    Print statement to keep track of things (slows things down!)
	    write(*,*) 'Iteration:', i_inc, ' | Alpha:', alpha, ' | Omega:', omega/(pi), ' | Inc:', i_p/pi

	    w_peri=0.0
	    ecc=0.0 
	    phi_0=2.0d0*pi*ran1(idum)    
	    mag_max=0.0d0
	    t_max=0.0d0
	    phi_mag_max=0.0d0   
	    nrnew=3 

c	    We express distances in terms of the Einstein angle of VB 10,
c	    which we have taken to be 10 mas.
	    mag_1=1.0d0  
	    mag_2=1.0d0
	    dev1=0.0
	    dev2=0.0 
	    mag_dif=0.0d0 
	    t_dif=0.0d0 
	    mag_dif_max=0.0d0 
	    phi_dif_mag_max=0.0d0

c	    t_early and t_late are latest and earliest times of peaks. counters count the numebr of peaks and valleys >0.01 deviation
	    t_early=1000000000.0d0
	    t_late=0.0d0
	    counter_peak=0
	    counter_vall=0

c	    arrays of integers reset to 0
	    Do 101 iarr=1,7
	      up(iarr)=0
	      down(iarr)=0
	      t_pre(iarr)=0
	      t_post(iarr)=0
	      magup(iarr)=0
	      magdown(iarr)=0
101	    enddo
	    caust=0 
	    val0001=0
	    val0001=0

c	    now do light curves 
	    sample=1
	    Do 102 it=1,tmax
	      caust1=0  
	      time=real(it)
              phi=phi_0-(time/P)*(days(it)-1) 	
	      xcm=xcms(it)
	      ycm=beta

c	coordinates of two binary objects. 
              x1=xcm-a1*((dcos(omega)*dcos(phi))-
     *	        (dsin(omega)*dsin(phi)*dcos(i_p)))
	      y1=ycm-a1*((dsin(omega)*dcos(phi))+
     *		(dcos(omega)*dsin(phi)*dcos(i_p)))
	      x2=xcm+a2*((dcos(omega)*dcos(phi))- 
     *		(dsin(omega)*dsin(phi)*dcos(i_p)))
	      y2=ycm+a2*((dsin(omega)*dcos(phi))+
     *		(dcos(omega)*dsin(phi)*dcos(i_p)))
	      u2=dsqrt(x2**2+y2**2)
	      mag_bin=mag(0.0d0,0.0d0,x1,y1,x2,y2,m1/mtot,m2/mtot,nrnew)
	      if(i_inc.eq.1) write(10,*) time, mag_bin-1.0, x1, y1, x2, y2 
	      dev_bin=mag_bin-mag_pt(it)

c	      (dev_bin is deviation from point lens)
	      if (nrnew.eq.5) caust=1
	      if (nrnew.eq.5) caust1=1

c-----------------------------------------------------------------------------------------
c	    JM: this next section might be a little tricky to follow- it identifies peaks in the magnification, 
c	    but not particularly well.
c
c           is there a peak in the deviation? this section finds out and bins integers according to locations/mags/devs
	      If(it.ge.3*sample)Then
		if(dev_abs.ge.0.001)then
		  If((dev1.lt.dev2).and.(dev1.lt.dev_bin)) val0001=val0001+1
		  If((dev1.gt.dev2).and.(dev1.gt.dev_bin)) peak0001=peak0001+1
		endif	
		dev1=mag_1-mag_pt(it-sample)
	        dev2=mag_2-mag_pt(it-2*sample)
	        dev_abs=dabs(dev1)


	        if(dev_abs.ge.0.01) then 		
c		dev1 and dev2 are the deviation corresponding to mag1 and mag2 as before. so we study dev1 and see if it is a maxima. 

c		is there a valley?
c		if there is a valley with deviation >0.01 then we bin according to characteristics, and count the valley

	          If((dev1.lt.dev2).and.(dev1.lt.dev_bin))Then 
	 	    counter_vall=counter_vall+1
		    t_dif=days(it-sample)
		    phi_dif_mag_max=phi
		    tdev=t_dif-t_pt_max	
		   write(*,*) 'PEAK!'	  
		    write(40+imass,900) -1, ialpha, alpha, p_day, phi_0, omega, i_p, nrnew, dev1, mag_1, t_dif, tdev             	            
		  endif


c		is there a peak?
c		if there is a peak with deviation >0.01 then we bin according to characteristics, and count the peak

		  If((dev1.gt.dev2).and.(dev1.gt.dev_bin))Then 
	 	    counter_peak=counter_peak+1
		    t_dif=days(it-sample)
		    phi_dif_mag_max=phi
c		  early and late times of peak, tdev is time deviation from pt lens max
		    if (days(it).gt.t_late) t_late=days(it)
		    if (days(it).lt.t_early) t_early=days(it)
	            tdev=t_dif-t_pt_max		
		    write(*,*) 'VALLEY!'  
		    write(40+imass,900) 1, ialpha, alpha, p_day, phi_0, omega, i_p, nrnew, dev1, mag_1, t_dif, tdev 	          	    	  	
		  endif

	        endif

	      endif

c---------------------------------------------------------------------------------------------

c	      now set mag2,1 to next ones along, and test if the maximum magnification is exceeded

	      mag_2=mag_1
	      mag_1=mag_bin
	      If(mag_bin.gt.mag_max)then
	        mag_max=mag_bin
	        t_max=days(it) 
	        phi_mag_max=phi
	      End If 
	      u1=dsqrt(x1**2+y1**2) 
	      u2=dsqrt(x2**2+y2**2)/dsqrt(m1/mtot)  
	      mag_pt1=(u1**2+2.d0)/(u1*sqrt(u1**2+4.d0))  
	      mag_pt2=(u2**2+2.d0)/(u2*sqrt(u2**2+4.d0))
	  
c	    end of one light curve:  
102	    End Do 

	    if (t_early.eq.1000000000.0) t_early=t_pt_max
	    if (t_late.eq.0.0) t_late=t_pt_max
c------------------------------------------------------------------------------------------------------------	
 	
c	finished looping over the light curve, so write arrays of integers for that light curve to file.
c	write all the info to a file for that sepcific planet and path.
c	format files: 900 is for peaks, 901 is for light curve.
c	alpha, p_day, phi_0, omega, i_p: alpha, period in days, initial phase, omega, inclination
c	caust: 1 if caustic crossing, 0 without. mag_max is mac magnification. val0001, peak0001 number of peaks and valleys 
	  
	    write(30+imass,901) alpha, p_day, phi_0, omega, i_p, mag_max, caust, val0001, peak0001, t_early, t_late, counter_peak, counter_vall

c	format staements for two filenames.
900	format(I4, I10, F14.6, F14.6, F14.6, F14.6, F14.6, I4, F14.6, F14.6, F14.6, F14.6)
901	format(F14.6, F14.6, F14.6, F14.6, F14.6, F14.6, I4, I4, I4, F14.6, F14.6, I4, I4) 

c	end of looping over orbital elements and alpha: 
100	End Do  

c	end of looping over planet masses, close path file...
	close(30+imass)
	close(40+imass)

c	call infowrite_end subroutine to write out useful stuff on screen
	call infowrite_end(lc_filename, peak_filename)


	End	  
c	END OF ACTUAL MOVING CODE. 





c	FUNCTIONS AND SUBROUTINES BELOW

c----------------------------------------------------------------------------------------------------------
c	define functions. ran1 is a random number generator. tmf is stellar mass function. lum is luminosity from core mass.

        function tmf(ran)
	Implicit none
	Real tmf,ran  
c       This function chooses a stellar mass from the Miller-Scalo
c       initial mass function.
        tmf=0.19*ran/((1.-ran)**.75+.032*(1.-ran)**.25)
        end

        function lum(em,cm)
         Implicit none
         Real lum,cm,em
c       lum is the luminosity as calculated by using the core mass
         lum=(2.d5*cm**6)/(1.d0 + 2.5d0*cm**4 + 3.d0*cm**5) +
     *             (2.2d0*em**3.2d0)
        end

c	JM/ I don't know what this next function does
        function cof(ms)
        Implicit none
        Real cof,ms,co0
        cof=10.**(-0.22+0.36*(log10(ms+1.))**2.5)
        end

c	RNG- requires initialising at start of code (seeding)
        function ran1(idum)
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        REAL ran1,AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        DATA iv /NTAB*0/, iy /0/
        if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
11        continue
          iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran1=min(AM*iy,RNMX)
        return
        end

c	function that calculates magnification of single lens

        function amp(u)
         implicit none
         real*8 amp,u
         amp=(u**2+2.d0)/(u*sqrt(u**2+4.d0))
        end


c------------------------------------------------------------------------------------------------------------
c       calculate the total magnification for a binary lens
c       xs, ys:, source position
c       m1, x1, y1: mass and position for the first lens
c       m2, x2, y2: mass and position for the second lens
c------------------------------------------------------------------------------------------------------------

c	MAIN MAGNIFICATION FUNCTION: solves 5th order polynomial, computes magnification for specific arrangement of lens and source.

        function mag(xs, ys, x1, y1, x2, y2, m1, m2,nr)

        implicit none
        real*8 mag
        real*8 xs, ys
        real*8 x1, y1, x2, y2
        real*8 m1, m2
        complex*16 z2, z2c
        complex*16 zs, zsc
        complex*16 z, zc
        complex*16 dzs
        integer i, m
        parameter (m=5)

c        common /nroots/ nr
        integer nr
        complex*16 b(m+1), roots(m), delta

c       eps solution error precision (should be larger for large mass ratio)
        real*8 eps
        parameter (eps=1.0d-5)

c       do a translation such that the first lens is at the origin
        zs = dcmplx(xs-x1, ys-y1)
        zsc = dconjg(zs)
        z2 = dcmplx(x2-x1, y2-y1)
        z2c = dconjg(z2)


c       polynomial coefficients generated by Mathematica

        b(6) = -zsc**2 + zsc*z2c
        b(5) = -(m1*zsc) - m2*zsc + zs*zsc**2 + 2*zsc**2*z2 + m2*z2c
     &     - zs*zsc*z2c -  2*zsc*z2*z2c
        b(4) = 2*m1*zs*zsc + 2*m2*zs*zsc + 2*m1*zsc*z2 - 2*zs*zsc**2*z2
     &     - zsc**2*z2**2 - m1*zs*z2c - m2*zs*z2c - m2*z2*z2c
     &     + 2*zs*zsc*z2*z2c + zsc*z2**2*z2c
        b(3) =  m1**2*zs + 2*m1*m2*zs + m2**2*zs - m1*m2*z2 - m2**2*z2
     &     - 4*m1*zs*zsc*z2 - 2*m2*zs*zsc*z2 - m1*zsc*z2**2
     &     + m2*zsc*z2**2 + zs*zsc**2*z2**2 + 2*m1*zs*z2*z2c
     &     + m2*zs*z2*z2c - zs*zsc*z2**2*z2c
        b(2) = -2*m1**2*zs*z2 - 2*m1*m2*zs*z2 + m1*m2*z2**2
     &     + 2*m1*zs*zsc*z2**2 - m1*zs*z2**2*z2c
        b(1) = m1**2*zs*z2**2


c       solve the polynomial equation, see Numerical Recipies

        call zroots(b, m, roots, .true.)


c       check whether the solution is a real solution of the lens equation

        mag = 0.0d0
        nr = 0
        do i=1, m
           z = roots(i)
           zc = dconjg(z)
           delta = zs-z+m1/zc + m2/(zc - z2c)

c          real solution, calculate the magnification

           if (cdabs(delta) .lt. eps) then
                nr = nr + 1
                dzs = m1/z**2+m2/(z-z2)**2
                mag = mag + 1.0/cdabs(1.0 - dzs*dconjg(dzs))
           endif
        enddo

        return
        end

c-------------------------------------------------------------------------

c	zroots subroutine, called in mag

      SUBROUTINE zroots(a,m,roots,polish)
      INTEGER m,MAXM
      REAL*8 EPS
      DOUBLE COMPLEX a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (EPS=1.d-10,MAXM=101)
CU    USES laguer
      INTEGER i,j,jj,its
      DOUBLE COMPLEX ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=dcmplx(0.d0,0.d0)
        call laguer(ad,j,x,its)
        if(dabs(dimag(x)).le.2.d0*EPS**2*dabs(dreal(x)))
     &          x=dcmplx(dble(x),0.d0)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),its)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(dreal(roots(i)).le.dreal(x))goto 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue
      return
      END
c-------------------------------------------------------------------------

c	laguer subroutine, called in zroots

      SUBROUTINE laguer(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      DOUBLE COMPLEX a(m+1),x
      PARAMETER (EPSS=1.d-10,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      DOUBLE PRECISION abx,abp,abm,err,frac(MR)
      DOUBLE COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=cdabs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=cdabs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=cdabs(b)+abx*err
11      continue
        err=EPSS*err
        if(cdabs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.d0*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=cdabs(gp)
          abm=cdabs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.d0) then
            dx=m/gp
          else
            dx=cdexp(dcmplx(dlog(1.0d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END   


c----------------------------------------------------------------------------------------------------------

c	JM/ simple subroutines for writing out info on screen to improve above code look.

	subroutine infowrite()
	  implicit none
	  write(*,*) '###############################################################'
  	  write(*,*) '\n PROGRAM: Binary Lens Light Curves for General Orbital Elements\n'
  	  write(*,*) 'This code is intended to be an all encompassing fortran code for ' 
   	  write(*,*) 'elliptical and inclined orbits. We compute light curves for a large' 
	  write(*,*) 'number of different types of orbits, and the code identifies ' 
	  write(*,*) 'maxima in the deviation from the standard point lens light curve.\n' 

	  write(*,*) 'Step 1: Specify general parameters, e.g masses, lens distance, PM, beta'
	  write(*,*) 'Step 2: Iterate over randonly generated phase and orbital elements'
	  write(*,*) 'Step 3: Compute light curve along path, and identify peaks'
	  write(*,*) 'Step 4: Write peaks and light curves to file\n'
	  write(*,*) '###############################################################\n'
	return
	end

	subroutine infowrite_end(file1, file2)
	  implicit none
	  character(15) file1, file2
	  write(*,*) 'Light Curve computation complete. Your filenames the data is saved in are:'
	  write(*,*) file1, file2
	  write(*,*) 'Format for', file1,': alpha, p_day, phi_0, omega, i_p: alpha, period in days, initial phase, omega, inclination'
	  write(*,*) 'caust: 1 if caustic crossing, 0 without. mag_max is max magnification. val0001, peak0001 number of peaks and valleys '
	  write(*,*) 'Format for', file2,': 1/-1, ialpha, alpha, p_day, phi_0, omega, i_p, nrnew, dev1, mag_1, t_dif, tdev'
	  write(*,*) 'Records all peaks and their characteristics.'
	  write(*,*) '###############################################################\n'
	return
	end

	subroutine help()
	  write(*,*) 'This is some help'
	return
	end
c	/JM End of everything
c	all done.
c----------------------------------------------------------------------------------------------------------

   
