c	This program creates the 5x3 plot required
	Implicit none
	Real*8 pi,x,y,mtot,amp,u     
	Real*8 msun, c, G, pc, au
	Parameter(pi=3.14159265, G=6.67d-11, c=299792458.0d0, msun=1.98892d30, pc=3.086d16, au=1.5d11) 
	Integer i,j,k, it,idum,nr,nr0,ia,iq,tmax, t, nref, num, itmax, al, num5, num10,i_t, i_planet, n_run, i_run, ittt
	Integer num2, num4, num8, num16, num20, num25, num32, ilc, nrnew
	Real*8 a,mag,incx, incy, trans, yrs, xdat, ydat, distance, uu, xcm, ycm, A_b, A_s, dif, teta_e, psi_inc, psi_rot  
	Real*8 q, alpha,m1,m2,a1,a2, P, time, phi, aa, timereal, einstein, phi_0, time_max, peak_max, clock, u_pl, phi_max
	Real*8 x1,x2,y1,y2, s1, s2, dev, devref, u_cm, n_m, pro_ret, prob_tot, psi_incc, psi_rott, increment, u_ref
	Real*8 peak, prob, tpeak, peakdev, p_day, psi_y, psi_x, omega
	Real*8 xpeak,ypeak,xpeakplanet,ypeakplanet, n_mj, ph, alpha2, v, i_p, x_i, mu, beta, d_l, omega2, p_hour
        Real ran1
	Parameter(tmax=48000)
	real*8 xpos(tmax),ypos(tmax), times(tmax), probs(tmax), probs5(tmax), days(tmax)
	character(3) pathlabel
	character lab

c--------------------------------------------------------------------------------------------------------------------------
	
	write(*,*) 'PROGRAM: working out a light curve for an inclined orbit, esin paper 2012'


c	name your light curve file here
	open(unit=30,file='lcinc_opi4_ipi4_b0.5')

c	here you specify the values of
c	v velocity in m/s, the two masses m2 and m1 in solar masses, the distance d_l in pc to the lens
c	the two inclination angles, omega and i_p (i_p=0.5pi is edge on), and the value of beta, clock specifies anticlockwise or clockwise orbit
	m2=1.3
      m1=0.3
	alpha=0.35 
	d_l=50.0
	v=10000.0d0
	
c	These next two values determination the orientation of the system, i_p being inclination and omega being ascending node longitude
	i_p=0.25d0*pi
	omega=0.25d0*pi
	
c	beta is the angular distance of closest approach	
	beta=0.5d0
	
c	variable I used to use to vary clockwise or anti clockwise
c	now you can do this using the orbital orientation angles - if you flip an orbit round by pi radians
c	it will be going in the opposite direction	
	clock=1.0d0

c	the write statement looks like this:
c	write(30,*) time, A_b, A_b-A_s, xcm, ycm 
c	write the time in hours, binary lens mag, single lens mag, and center of mass coords

c------------------------------------------------------------------------------------------------------------------------------

c   total mass of system
	mtot = m1 + m2          
c 	mass ratio			
	q = m1 / m2
	
c	theta_e in arc seconds	
	teta_e = 2.0d0 * dsqrt((G*mtot*msun)/(d_l*pc*c**2)) * (180.0d0/pi) * 3600.0
	
c	proper motion in as / yr
	mu=((v/au)*(3600.0d0))/(d_l*teta_e)
	
c	einstein radius in 	
	einstein=teta_e*d_l
	
c	physical separation	
	aa=alpha*einstein	

c	orbital period, days		 
	P=365.25d0*dsqrt(aa**3/mtot)
	
c 	orbital period, hours	
	p_hour=P*24.0
	
c	initial phase, somewhat arbritrary
	phi_0=pi*0.5d0
	
c	initial position
	x_i=-( mu * 24000.0d0 ) + 0.000000001
	
c-------------------------------------------------------------------------------------------------------------------------
        a1=alpha * m2/mtot
        a2=alpha * m1/mtot
	itmax=1
	peakdev=0.0
	omega2=2.0d0*pi/p_hour
	Do 503 it=1,1	 
          a1=alpha*m2/mtot
          a2=alpha*m1/mtot
	  write(*,*) a1, a2
	  Do 504 i=1,tmax
	    time=real(i)
	    phi=phi_0-omega2*clock*(time) 
	    xcm=x_i+(mu*time)
	    ycm=beta
	    x1=xcm-a1*((dcos(omega)*dcos(phi))-
     *	      (dsin(omega)*dsin(phi)*dcos(i_p)))
	    y1=ycm-a1*((dsin(omega)*dcos(phi))+
     *	      (dcos(omega)*dsin(phi)*dcos(i_p)))
	    x2=xcm+a2*((dcos(omega)*dcos(phi))- 
     *	      (dsin(omega)*dsin(phi)*dcos(i_p)))
	    y2=ycm+a2*((dsin(omega)*dcos(phi))+
     *	      (dcos(omega)*dsin(phi)*dcos(i_p)))
	    if (i.eq.24000) write(*,*) x, y, x1, y1, x2, y2, omega2, time
	    u_cm=dsqrt(xcm**2+ycm**2)
	    A_s=(u_cm**2+2.d0)/(u_cm*sqrt(u_cm**2+4.d0))
	    x=0.0d0
	    y=0.0d0 
	    A_b=mag(0.0d0,0.0d0,x1,y1,x2,y2,m1/mtot,m2/mtot,nrnew)
	    u_pl=dsqrt(x1**2+y1**2)
	    write(30,*) time, A_b, A_b-A_s, xcm, ycm        
504	  end do
503	End Do
	End   
c-------------------------------------------------------------------------------------------------------------------------
 

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
c       lum is the luminosity as calcbulated by using the core mass
         lum=(2.d5*cm**6)/(1.d0 + 2.5d0*cm**4 + 3.d0*cm**5) +
     *             (2.2d0*em**3.2d0)
        end

        function cof(ms)
        Implicit none
        Real cof,ms,co0
        cof=10.**(-0.22+0.36*(log10(ms+1.))**2.5)
        end

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
11      continue
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
      END

*******************************************************************************
c       calcbulate the total magnification for a binary lens
c       xs, ys:, source position
c       m1, x1, y1: mass and position for the first lens
c       m2, x2, y2: mass and position for the second lens
c
********************************************************************************





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

c
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
c
c       solve the polynomial equation, see Numerical Recipies

        call zroots(b, m, roots, .true.)
c
c       check whether the solution is a real solution of the lens equation

        mag = 0.0d0
        nr = 0
        do i=1, m
           z = roots(i)
           zc = dconjg(z)
c           write(6, *) i, roots(i)+dcmplx(x1, y1)
           delta = zs-z+m1/zc + m2/(zc - z2c)
c
c          real solution, calcbulate the magnification
c          real solution, calcbulate the magnification

           if (cdabs(delta) .lt. eps) then
                nr = nr + 1
                dzs = m1/z**2+m2/(z-z2)**2
                mag = mag + 1.0/cdabs(1.0 - dzs*dconjg(dzs))
           endif
        enddo

        return
        end

        function amp(u)
         implicit none
         real*8 amp,u
         amp=(u**2+2.d0)/(u*sqrt(u**2+4.d0))
        end

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

   

