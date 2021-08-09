!              The main driver for the solution of the kinetic 
!              equations describing particle injection and energy losses,
!              including photopair (Bethe-Heitler) and photomeson production from SOPHIA.
!              The NAG library routine D02EJF is used (https://www.nag.com/numeric/fl/nagdoc_fl24/pdf/d02/d02ejf.pdf)

!              SUMMARY OF PHYSICAL PROCESSES
!              Electrons/positrons: 
!              external injection, physical escape, losses due synchrotron radiation, 
!              losses due to inverse Compton scattering, injection by photon-photon pair production,
!              injection from Bethe-Heitler pair production, injection from photomeson production.
!              Photons: 
!              production from synchrotron radiation, attenuation from synchrotron self-absorption, 
!              production from inverse Compton scattering, injection from photomeson production, 
!              attenuation from photon-photon pair production, physical escape.
!              Protons:
!              external injection, physical escape, losses due to synchrotron radiation 

!              There are 9 FLAGS in code.inp controling these processes:
!              isyn=0/1   switches off/on electron  synchrotron cooling.
!              iprsyn=0/1   switches off/on proton synchrotron cooling. 
!              icompt=0/1 switches off/on electron Compton cooling.
!              ikn=0/1    switches off/on electron Compton cooling in KN regime.
!              igg=0/1    switches off/on photon-photon pair production.
!              issa=0/1   switches off/on synchrotron self absorption.
!              iesc=0/1   switches off/on electron escape.
!              ipsc=0/1   switches off/on photon escape. 
!              ianni=0/1  switches off/on electron-positron annihilation. 

!              USEFUL DIMENSIONLESS PHYSICAL QUANTITIES.  
!              TB: magnetic compactness defined as (B**2/8π)/σΤ*R, with R the size of the blob.
!              Q: the ratio of the magnetic field to the critical value Bcr=4.413*10**13 G

!              DIMENSIONS OF CODE QUANTITIES.
!              Time is measured in units of the photon crossing time R/c.
!              Particle number densities are in units of σT*R.
!              Photon energies are in units of me*c**2.


!              PARTICLE ARRAYS.
!              YP(NP): array of dimension NP containing the natural log of the dimensionless 
!                      proton number density in equally spaced intervals of deltap=d(logγ).    
!              GP(NP): array of dimension NP containing the log10 of proton Lorentz factor γ. 
!              YE(NE): array of dimension NE containing the natural log of the dimensionless 
!                      electron/positron number density in equally spaced intervals of deltap=d(logγ).   
!              GE(NE): array of dimension NE containing the log10 of electron/positron Lorentz factor γ. 
!              YG(NG): array of dimension NG containing the natural log of the dimensionless 
!                      photon number density in equally spaced intervals of deltax=d(logx)=2*deltap.  
!              XG(NG): array of dimension NG containing the log10 of the dimensionless photon energy x.
!                      The lower and upper limits of the array are: X(1)=XMIN=Q*GP(1)**2 and X(NG)=XMAX=GP(NP).
!                      XSCH=X(NP)=Q*GP(NP)**2 is the maximum energy of synchrotron photons, so there
!                      is no synchrotron contribution from XSCH to XMAX. The number of photon extra bins 
!                      NADD=NG-NP is given by NADD=LOG(XMAX/XSCH)/DELTAX. NADD must be integer, so 
!                      care must be taken when choosing GP(1), GP(NP), and NP.

      implicit real*8 (a-h,o-z)

      integer n,iw,ifail
      integer iprext,ielext,ielextbr, iprextbr
      integer iphotext, iphotext2
      integer isyn,icompt,igg,issa,iesc,ipsc
      integer ianni
      integer tarray
      real*8 npdec
      character*1 relabs
	
!     PARAMETERS USED IN d02ejf ROUTINE 
      PARAMETER(npmax=400,ntotal=800,nwork=440350,relabs='M')
	
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
             
!     ARRAY DEFINITIONS
      DIMENSION tarray(3)
      DIMENSION y(ntotal),yp(npmax),ye(npmax),yg(npmax)
     $         ,gp(npmax),ge(npmax),x(npmax),work(nwork)
     $         ,ygamma(ntotal),ygbb(npmax) 

!     DECLARATION OF SUBROUTINES      
      external out
      external deriv
      external d02ejf,d02ejy,d02ejw

      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/flags/isyn,iprsyn,icompt,igg
     $  ,issa,iesc,ipsc,ianni,ikn 
      common/ygm/ygamma
      common/nonlin/factor
      common/tfl/xnorm
      common/htm/h,tend,iout  
      common/fqq/nsteps,nout,ireadin
      common/phexter/gext,ygext,iphotext  
      common/ext/x1,x2,xbr,beta1,beta2,extph0,iphotext2 
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp 
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,iprext,ipexp  
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr 
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr 
      common/ebb/ygbb,nbbmx	
      common/tvv/iprcl 
      common/pgp/radius,xl10min,bfield
      
      call cpu_time (t_start)
 
      call itime(tarray)
      tt1=tarray(1)*3600.+tarray(2)*60.+tarray(3)
 
      open (unit=13, file='code_hao.inp', status='old')
      open (unit=98, file='code_r4.dat', status='unknown')       
      open (unit=99, file='dum_r4.dat', status='unknown') 

!     READING THE INPUT FILE
      read (13,*) ireadin,npdec,nsteps,nout,tend
      read (13,*) gpexmx,ypmin,yemin,ygmin,tol
      read (13,*) slpinj,sleinj,slginj
      read (13,*) radius, bfield
      read (13,*) iprext,gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,ipexp
      read (13,*) ielext,geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ieexp
      read (13,*) iphotext,temperat,exlumth
      read (13,*) iphotext2,x1,x2,xbr,beta1,beta2,extph0
      read (13,*) ielextbr,geextbr,slelints1,slelints2,ae2	
      read (13,*) iprextbr,gpextbr,slprints1,slprints2,ap2
      read (13,*) isyn, iprsyn, issa 
      read (13,*) icompt, ikn, igg, ianni, iesc, ipsc
 
      
      if (ielextbr.gt.0) ielext=0      	
      if (iprextbr.gt.0) iprext=0	
 
	geextbr=10.**geextbr
	gpextbr=10.**gpextbr		
	gpextmn=10.**gpextmn
	gpextmx=10.**gpextmx
	geextmn=10.**geextmn
	geextmx=10.**geextmx
	exlumel=10.**exlumel
	exlumpr=10.**exlumpr

!       dimensionless B field	
	q=bfield/Bcr
!       magnetic energy density [erg/cm^3]	
	UB=Bfield**2./8./pi
!       magnetic compactness	
	tb=sthom*UB*radius/elrms
!       co-moving isotropic-equivalent electron injection luminosity [erg/s]
        xLeinj=exlumel*4.*pi*radius*elrms*c/sthom
!       co-moving isotropic-equivalent proton injection luminosity [erg/s]
        xLpinj=exlumpr*4.*pi*radius*prms*c/sthom
          
        write(6,*) 'B [G], UB [erg/cm^3], R [cm]'
	write (6,1000) Bfield,UB,Radius
	write(6,*) 'Le,inj [erg/s], Lp,inj [erg/s]'
        write(6,1000) xLeinj, xLpinj 
        
!       For blackbody	
	xbb=temperat*boltz/elrms
	yglbb=exlumth
	bblum=4*pi*radius*elrms*c/sthom
	Ubb=bblum/4./pi/radius**2/c
	Ubbrl=astbo*temperat**4

!       Normalization factors         
        factor=tb*q**(-5./3.)/1836.1
	xnorm=tb*q**(-5./3.)
	
	
      if (ireadin.eq.0) then
      write(98,1002)ireadin,npdec,nsteps,tend
      write(98,1005)gpexmn,gpexmx,xexmn,tol 
      write(98,1005)slpinj,sleinj,slginj
      write(98,1000)bplus,tb,q,radius
      write(98,1000)ypmin,yemin,ygmin
      write(98,1001)iphotext,temperat,exlumth
      write(98,1001)iphotext2,x1,x2,xbr,beta1,beta2,extph0
      write(98,1003)isyn, iprsyn, issa
      write(98,1003)icompt, ikn, igg, ianni, iesc, ipsc      
      end if



!     LIMITS FOR THE PARTICLE ENERGY ARRAYS -- NOT TO BE CHANGED
        gpexmn=1./(1.*npdec)
	iprcl=0

	gpmin=10.**gpexmn
	gpmax=10.**gpexmx	
	gemin=gpmin
	gemax=100.*gpmax

	
!       total number of bins for the protons/electron arrays
	np=int(npdec*log10(1.01*gpmax))
	ne=int(npdec*log10(1.01*gemax))
	
!       logarithmic step for particle and photon arrays
        deltap=log(gpmax/gpmin)/real(np-1)
	deltax=2.*deltap
	
!       change in the binning of the photon array for high B.
        if(bfield.GT.1.E2)then    
          abin = int(log10(bfield/1.E2))+1
          xmin=q*gemin**2./(10**abin)
          else
	  xmin=q*gemin**2.
        end if
	xlmin=log(xmin)
	xl10min=log10(xmin)
	
!       determine the dimension of the photon array NG
	do n=1,3*ne
	 x(n)=exp(xlmin+(n-1)*deltax)
	 if (x(n).gt.gemax) then
	 ng=n-1
	 goto 300
	 endif 
        enddo
300     continue
	xmax=x(ng)
!         ntot=np+ne+ng+1+ne+np	!nt!N
        ntot=np+ne+ng+1

      write (6,*) 'Np=',np
      write (6,*) 'Ne=',ne
      write (6,*) 'Ng=',ng
      write (6,*) 'Ntot=',ntot

      write(6,7000)gpmin,gpmax,gemax
      write(6,7001)xmin,xmax

7000  format(1x,'min/max limits for protons and electrons:',3(g12.5,1x))
7001  format(1x,'min/max limits for photons:',2(g12.5,1x))


!     INITIALIZATION OF THE PARTICLE ARRAYS
      yp(1)=ypmin*log(10.)
      y(1)=yp(1)
      ye(ne)=yemin*log(10.)
      y(np+ne)=ye(ne)
      yg(ng)=ygmin*log(10.)
      y(np+ne+ng)=yg(ng)
      gp(np)=gpmax
      ge(ne)=gemax


!     background protons/neutrons
      do n=1,np
         gp(n)=gpmin*(gpmax/gpmin)**(real(n-1)/real(np-1))
         yp(n)=(-slpinj)*log(gp(n)/gp(1))+ypmin*log(10.)-20.         
         y(n)=yp(n)         
         ygamma(n)=gp(n)
      enddo

!     background electrons/muons/pions/kaons     
      do n=ne-1,1,-1
         ge(n)=gemin*(gemax/gemin)**(real(n-1)/real(ne-1))
         ye(n)=(ne-n)*sleinj*deltap + ye(ne)
         y(n+np)=ye(n)        
      enddo 
         
	 do n=1,ne
	 ygamma(np+n)=ge(n)
	 enddo 
      
!     background photons 
!     ygbb: the external photon field dimensionless number density
!     yglbb: external photon compactness for BB field (read from input file)
!     extph0: external photon compactness for BPL field (read from input file)
!     we set initially ygbb<<yg

!!!   external photons with blackbody (BB) energy spectrum 
      if(iphotext.eq.1)then 
      sum=0.
      do n=1,ng
         yg(n)=-slginj*log(x(n)/x(ng)) + yg(ng)
	 ygbb(n)=yg(n)-20.
            if (x(n).gt.30.*xbb) goto 301  
	bbnorm=45.*yglbb/pi**4./xbb**4./xnorm 
 	ygbb(n)=log(bbnorm*x(n)**2./(exp(x(n)/xbb)-1.))
	sum=sum+deltax*x(n)**2.*exp(ygbb(n))
301     continue  
	y(np+ne+n)=log(exp(yg(n))+exp(ygbb(n))) 
        enddo
	endif
	
!!!     external photons with broken power law (BPL) spectrum 
        if(iphotext2.eq.1) then
        sum=0.
        do  n=1,ng
        yg(n)=-slginj*log(x(n)/x(ng)) + yg(ng)
        ygbb(n)=yg(n)-20.
           if(x(n).ge.x1.and.x(n).le.xbr)then
           ygbb(n)=log(x(n)**(-beta1)/xnorm)
           elseif(x(n).gt.xbr.and.x(n).le.x2)then
           ygbb(n)=log(xbr**(beta2-beta1)*x(n)**(-beta2)/xnorm)
           end if
        sum=sum+deltax*x(n)**2.*exp(ygbb(n))
        enddo
!      find correct normalization for BPL 
       bbnorm=3*extph0/(xnorm*sum) !! correction from code comparison project -- 27/04/2020 !!
       sum=0.
!      add to the photon bg the BPL field with correct normalization       
       do n=1,ng
       yg(n)=-slginj*log(x(n)/x(ng)) + yg(ng)
       ygbb(n)=yg(n)-20.
       if (x(n).gt.3.d0*x2) then
       nbbmx=n-1
       endif
           if(x(n).ge.x1.and.x(n).le.xbr)then
	   ygbb(n)=log(bbnorm*x(n)**(-beta1)/xnorm)
	   elseif(x(n).gt.xbr.and.x(n).le.x2)then
	   ygbb(n)=log(bbnorm*xbr**(beta2-beta1)*x(n)**(-beta2)/xnorm)
	   endif
       sum=sum+deltax*x(n)**2.*exp(ygbb(n))	   
       y(np+ne+n)=log(exp(yg(n))+ exp(ygbb(n)))
       enddo
       endif

!   no external photons	
	if(iphotext.eq.0.and.iphotext2.eq.0)then
	do n=1,ng
	yg(n)=-slginj*log(x(n)/x(ng)) + yg(ng)
	y(np+ne+n)=yg(n)
	enddo
	endif
		
	do  n=1,ng
        ygamma(np+ne+n)=x(n)
        enddo 
 
        
	do n=1,ng
	if (iphotext.eq.1)then
	   if (x(n).gt.30.*xbb) then
	   nbbmx=n-1
	   goto 302
	   endif
	elseif (iphotext2.eq.1)then
	   if (x(n).gt.3.d0*x2) then
	   nbbmx=n-1
	   goto 302
	   endif
	endif    
        enddo
302     continue  

! the last point is reserved for the cool electrons (γ<=1)
	ynecool=ye(1)  
	y(np+ne+ng)=-200. 
	     
 
! 	t=0.
! 	iout=nsteps
	
	t=0.
	iout=nsteps

!      In case the program crashes:
!             Pick up where we left off...

         if(ireadin.eq.1)then
         read(99,*) iout,t
         do ijk=1,ntot
            read(99,*)y(ijk)
	    if (ijk.le.np) y(ijk)=y(ijk)-0.0*log(10.)
         enddo
         end if
         
	h=(tend-t)/real(iout)


!       CALL TO MAIN ROUTINE FOR SOLVING THE SYSTEM OF EQUATIONS	
	ifail=0
	iw=nwork	
      call  d02ejf(t,tend,ntot,y,deriv,d02ejy,tol,relabs,out,
     $             d02ejw,work,iw,ifail)
         
 
         if(tol.le.0)then
                      print*,' warning, no change in solution'
                      tol=abs(tol)
         end if
         
         if(ifail.ne.0) goto 303

         print*,' integration complete to t=',tend
         call cpu_time(t_stop)
         write(*,*) 'Elapsed CPU time [s] = ', t_stop-t_start
         call itime(tarray)
         tt2=tarray(1)*3600.+tarray(2)*60.+tarray(3)
         write(*,*) 'Elapsed wallclock time [s] = ', tt2-tt1
      stop

!                   diagnostics for failed D02... call
 
303   continue

      write(6,6200)tol,tend,t
      write(6,6201)ntot
      
      do n=1,ntot
         write(6,6202)yg(n),y(n)
      enddo

1000  format(1x,6(1pe12.4,1x))
1005  format(3x,4(1pd12.4),2x,i10)
1006  format(2x,i10,2x,1pd12.4)
1002  format(2x,i10,2x,d12.4,2x,i10,2x,d12.4)
1001  format(1x,i10,6(1pd12.4))
1003  format(2x,13(i10))
1004  format(2x,2(i10),2x,2(1pd12.4))
1100  format(1x,g12.4,i5)
1200  format(1x,g12.4,i5,1x,1pe12.4,2x,e12.4)
6200  format(1x,'tol=',g12.4,' tend=',g12.4,' t=',g12.4)
6201  format(1x,'ntot=',i5)
6202  format(1x,e16.8,3x,e16.8)

!!!   END OF MAIN PROGRAM !!!      
      END
      
      
!!! SUBROUTINES !!! 
************************************************************************

      subroutine out(t,y)
!     OUTPUTS THE RESULTS      
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
     
      DIMENSION y(ntotal)
      DIMENSION yp(npmax),ye(npmax),yg(npmax),zd(npmax)
     $         ,gp(npmax),x(npmax),ygamma(ntotal),ge(npmax)
     $	       ,ygbb(npmax)

      common/param/tb,q,gpmax,
     $  deltap,deltax,np,ne,ntot,ng
      common/flags/isyn,iprsyn,icompt,igg
     $  ,issa,iesc,ipsc,ianni,ikn 
      common/ygm/ygamma 
      common/htm/h,tend,iout  
      common/tfl/xnorm
      common/fqq/nsteps,nout,ireadin
      common/ebb/ygbb,nbbmx
      common/pyy/tauth 
      common/xuu/xen,gmax
      common/pgp/radius,xl10min,bfield
 

      do n=1,np
	 gp(n)=ygamma(n) 
         yp(n)=y(n)
      enddo

      do n=1,ne
	 ge(n)=ygamma(np+n) 
         ye(n)=y(np+n)
      enddo

      do n=1,ng
        x(n)=ygamma(np+ne+n) 
	if (x(n).lt.0.1) zd(n)=1.
	if (x(n).ge.0.1.and.x(n).le.1.) zd(n)=(1.-x(n))/.9
	if (x(n).gt.1.) zd(n)=0.
        yg(n)=y(np+ne+n)-log(1.+tauth*zd(n)/3.)
      enddo
	

!  the one before last bin is reserved for the cooled electrons
        ynecool=ye(1) 
!  measure of the Thomson optical depth on cold electrons	
	tauth=exp(ynecool) 
	
! writes solution in the dum_r4.dat file. The code picks up from this point.
         rewind 99
         write(99,*)iout+1,t
         do  ijk=1,ntot
            write(99,*)y(ijk)
         enddo

!  sump is a measure of the energy density in protons
!  sumnp is a measure of the density in protons
!  sume is a measure of the energy density in electrons
!  send is a measure of the density in electrons
!  sumx is a measure of the energy density in photons
!  sxnd is is a measure of the density in photons

	   sump=0.
	   sumnp=0.
	   sume=0.
	   send=0.
           sumx=0.
	   sxnd=0. 
	   
	do k=1,np
	   if (k.eq.1.or.k.eq.np) then
	      fact=.5
		else 
	      fact=1.
	   end if	  
	   sump=sump+fact*deltap*gp(k)**2.*exp(yp(k))
           sumnp=sumnp+fact*deltap*gp(k)*exp(yp(k))
        enddo

        do k=1,ne
            if (k.eq.1.or.k.eq.ne) then
	      fact=.5
		else 
	      fact=1.
	   end if 
	   sume=sume+fact*deltap*ge(k)**2.*exp(ye(k))
	   send=send+deltap*ge(k)*exp(ye(k))
        enddo 
        
	 do k=1,ng
	   if (k.eq.1.or.k.eq.ng) then
	      fact=.5
		else 
	      fact=1.
	   end if 
	   sumx=sumx+fact*deltax*x(k)**2.*exp(yg(k))
	   sxnd=sxnd+fact*deltax*x(k)*exp(yg(k))
        enddo 
           xen=sumx*xnorm  

      write (6,1000) t,sump,xen,tauth
      write (90,1000) t,sump,xen,tauth 
      write (98,1000) t,sump,sume,sumnp,send,sxnd
      write (72,1000) t,log10(xen)
      write (73,1000) t,log10(ratmpe*sump) 
      
      write (81,1000) t,sump,xen,tauth	
	do n=1,ng	
	if(exp(yg(n))-exp(ygbb(n)).gt.0.)then
        write (71,1000) log10(x(n)),
     $  log10(xnorm*(exp(yg(n))-exp(ygbb(n)))*x(n)**2),
     $  log10(xnorm*exp(ygbb(n))*x(n)**2)
	else
    	write (71,1000) log10(x(n)),log10(xnorm*(exp(yg(n)))*x(n)**2),
     $	log10(xnorm*exp(ygbb(n))*x(n)**2)
	endif
	write (81,1000) log10(x(n)),log10(xnorm*exp(yg(n))*x(n)**2)	
        enddo

      write (88,1000) t,sump,xen,tauth
	do n=1,np
        write (88,1000) log10(gp(n)), log10(gp(n)**2*exp(yp(n)))
        enddo 
        
      write (89,1000) t,sump,xen,tauth
	do n=1,ne
        write (89,1000) log10(ge(n)), log10(ge(n)**2*exp(ye(n)))
        enddo  
	
      
      
	do n=1,ntot
	if (n.lt.ntot) then
	  if (n.le.np+ne) then
	    tau=0.
	    else
	    tau=tauth*zd(n-np-ne)/3.
	  end if
	endif
        enddo

	t=tend-real(iout)*h
	iout=iout-1
	
1000  format(1x,6(1pe12.4,1x))
	return
	end
	
**************************************************************************

      subroutine deriv(time,y,yprime)
!! COMPUTES TIME DERIVATIVES OF ALL PHYSICAL PROCESSES
!! subroutine called by d02ejf

      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer j1,j2,k1
      integer tarray
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,afine=7.3e-3,csrt=.2,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 

      DIMENSION y(ntotal),yprime(ntotal),ygamma(ntotal)
      DIMENSION yp(npmax),ye(npmax),yg(npmax),ge(npmax)
     $         ,gp(npmax),x(npmax),ygbb(npmax),yinit(npmax) 
     $         ,gcsth(npmax),gcskn(npmax),gdens(npmax)
     $         ,cskn(npmax),syncel(npmax),syncph(npmax)
     $         ,csth(npmax),denkn(npmax),zd(npmax),relann(npmax)
     $         ,sst(npmax),extelin(npmax),extprinxp(npmax)
     $         ,syncprph(npmax),syncpr(npmax)
      DIMENSION  yggal(npmax),yggcral(npmax)
     $   ,anng(npmax),annel(npmax)
     $   ,ggabsr(npmax),gginje(npmax)    
   
      DIMENSION tarray(3)
      
      common/param/tb,q,gpmax,
     $ deltap,deltax,np,ne,ntot,ng
      common/flags/isyn,iprsyn,icompt,igg
     $  ,issa,iesc,ipsc,ianni,ikn 
      common/ygm/ygamma 
      common/tfl/xnorm
      common/nonlin/factor      
      common/phexter/gext,ygext,iphotext 
      common/ext/x1,x2,xbr,beta1,beta2,extph0,iphotext2 
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp  
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,iprext,ipexp  
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr  
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr       
      common/ebb/ygbb,nbbmx 
      common/pyy/tauth 
      common/xuu/xen,gmax
      common/tvv/iprcl  
      common/fqq/nsteps,nout,ireadin
      common/pgp/radius,xl10min,bfield
 
	iprcl=iprcl+1 
      
      do n=1,np
	 gp(n)=ygamma(n) 
         yp(n)=y(n)
      enddo 

      do n=1,ne
	 ge(n)=ygamma(np+n) 
         ye(n)=y(np+n)
      enddo

	gemin=ge(1)
	gemax=ge(ne)

      do n=1,ng
	 x(n)=ygamma(np+ne+n) 
         yg(n)=y(np+ne+n)
	if (x(n).lt..1) zd(n)=1.
	if (x(n).ge..1.and.x(n).le.1.) zd(n)=(1.-x(n))/.9
 	if (x(n).gt.1.) zd(n)=0.
        enddo
	zd(ng)=0.
 
!  the one before last bin is reserved for the cooled electrons
        ynecool=ye(1) 
!  measure of the Thomson optical depth on cold electrons	
	tauth=exp(ynecool) 
	if (iprcl.eq.1) then
	n=np+ne+ng+1+ne+np	!N
	yinit(n)=y(n)
	endif

!!! RATES OF PHYSICAL PROCESSES !!!

!! PROTON INJECTION 
        call prinj(gp,yp,np,deltap,extprinxp,sumprinj)
!! ELECTRON INJECTION
        call elinj(ge,ye,ne,deltap,extelin,sumelinj)        
!! SYNCHROTRON RADIATION
        if(isyn.eq.1) then
! i) electron loss   
        call synecool(ge,ye,syncel,sumsyncel)
! ii) photons from electron synchrotron
        call synephot(ge,ye,x,yg,syncph,sumsyncph)
        endif 
        if(iprsyn.eq.1) then
! iii) proton losses
        call synpcool(gp,yp,syncpr,sumsyncpr)
! iv) photons from proton synchrotron        
        call synpphot(gp,yp,x,yg,syncprph,sumsyncprph)
        endif 
!! SYNCHROTRON SELF ABSORPTION
        if(issa.eq.1) then
! photon loss (attenuation)
        call ssaphot(ge,ye,x,yg,sst,sumphssa)
        else 
        do n=1,ng
        sst(n) = 0
        enddo 
        endif 
!! INVERSE COMPTON SCATTERING (ICS)        
        if (icompt.eq.1) then
! electron loss
        call compe_sim(x,ge,ye,yg,ikn,csth,cskn,sumcsth,
     $ sumcskn,taukn) 
! photons from electron ICS using Bloumenthal-Gould emissivity
        call compph_bg(x,ge,yg,ye,gcsth,sumgcsth)
        else 
        
        do n=1,ne
        csth(n) = 0
        cskn(n) = 0
        enddo 
        
        do n=1,ng
        gcsth(n) = 0
        enddo 
        
        endif
!! PHOTON-PHOTON PAIR PRODUCTION
        if (igg.eq.1) then 
! photon attenuation
        call ggphot(x,yg,ggabsr,sumggabs)       
! electron/positron injection 
        call ggelec(x,yg,ge,ye,gginje,sumggabs,sumgginj)  
        endif 
        
!!! THE EQUATIONS !!!
!       derivatives computed at first bin

!       protons:
        yprime(1)=extprinxp(1)-iesc*bpresc+syncpr(1)*iprsyn
!       electrons:
        yprime(np+1)=extelin(1)-iesc*belesc+syncel(1)*isyn+
     $  (csth(1) - cskn(1))*icompt + igg*gginje(1) 
!       photons:
 	yprime(np+ne+1)=-ipsc*1./(1.+tauth*zd(1)/3.)+
     $	isyn*syncph(1) + iprsyn*syncprph(1)+ issa*sst(1) +
     $  icompt*gcsth(1) -igg*ggabsr(1)   
!       cold electrons: 
	pycsth=4./3.*xnorm*gdens(1)*(gp(1)**2.-1.)*exp(ye(1))/tauth
	pysyn= 4./3.*tb*(gp(1)**2.-1.)*exp(ye(1))/tauth
        yprime(np+ne+ng+1)= pycsth+pysyn+taukn
     $    -3.*tauth/32.-belesc	     
 !       derivatives computed at last bin

!       protons: 
	yprime(np)=-.1
!       yprime(np)=extprinxp(np)-bpresc+syncpr(np)*isyn	
        if (y(np).lt.-150.) yprime(np)=0.
!       electrons:
	yprime(np+ne)=-.1
!       yprime(np+1)=extelin(ne)-belesc+syncel(ne)*isyn	
!        photons:
        yprime(np+ne+ng)=-ipsc*1./(1+tauth*zd(ng)/3.)+   
     $	isyn*syncph(ng) +iprsyn*syncprph(ng) +issa*sst(ng) +
     $  icompt*gcsth(ng) -igg*ggabsr(ng) 

***************************

!       derivatives computed at all other bins 

	sumprot=0.
	
!! PROTONS 
	do n=2,np-1!    
! protons	
	yprime(n)=extprinxp(n)-iesc*bpresc+iprsyn*syncpr(n)

        sumprot=sumprot+deltap*gp(n)**2.*exp(yp(n))*yprime(n)

        enddo 
    
        sumelec=0.
     
!! electrons       
        do n=2, ne-1 
        yprime(np+n)=extelin(n)-iesc*belesc+
     $  isyn*syncel(n)+
     $  icompt*(csth(n) - cskn(n)) +
     $  igg*gginje(n) 
	sumelec=sumelec+deltap*ge(n)**2.*exp(ye(n))*yprime(np+n)
        enddo 
  	
************************
!   PHOTONS
	
	sumphot=0.

	do n=1,ng-1
   
	yprime(np+ne+n)=-ipsc*1./(1.+tauth*zd(n)/3.) +
     $	isyn*syncph(n) +iprsyn*syncprph(n) +issa*sst(n) +
     $  icompt*gcsth(n) -igg*ggabsr(n)
  
        sumphot=sumphot+deltax*x(n)**2.*xnorm*exp(yg(n))*yprime(np+ne+n)

! treatment of external photon field     
        if ((iphotext2.eq.1.and.x(n).le.x2.and.x(n).ge.x1).
     $  or.(n.le.nbbmx.and.iphotext.eq.1)) then
	if (ygbb(n).gt.yg(n)) then
           yprime(np+ne+n)=-1./(1.+tauth*zd(n)/3.) + exp(ygbb(n)-yg(n))
	else
           yprime(np+ne+n)=yprime(np+ne+n)+exp(ygbb(n)-yg(n))
        endif
        endif 
        enddo
        

	yprime(np+ne+ng)=0.
	
! compute photon number density and photon energy density (dimensionless)	
        sumx=0.
	sumtx=0.
	 do k=1,ng
	   if (k.eq.1.or.k.eq.ng) then
	      fact=.5
		else 
	      fact=1.
	   end if
	   sumx=sumx+fact*deltax*x(k)**2.*exp(yg(k))	
	   tau=tauth*zd(k)/3.
	   sumtx=sumtx+fact*deltax*x(k)**2.*exp(yg(k))/(1.+tau)	
        enddo 
	xden=xnorm*sumx !-bblnum
	xdten=xnorm*sumtx !-bblnum/(1.+tauth/3.) 

    
      if(ireadin.eq.0.and.iprcl.eq.1)then ! BEGIN warnings for initial background choice  
        ypremx=abs(yprime(np+1))
        yprpmx=abs(yprime(1))
        yprgmx=abs(yprime(np+ne+1))

	nemx=1
	npmx=1
	ngmx=1

	do n=1, np
	if(yprpmx.lt.abs(yprime(n)))then
	yprpmx=yprime(n)
	npmx=n
	endif
	enddo
			
	do n=1, ne
	if(ypremx.lt.abs(yprime(np+n)))then
	ypremx=yprime(np+n)
	nemx=n
	endif
	enddo
        
        do n=1, ng
	if(yprgmx.lt.abs(yprime(np+ne+n)))then
	yprgmx=yprime(np+ne+n)
	ngmx=n
	endif
	enddo
	
        if(yprpmx.lt.1d10) then
        write(6,*) 'WARNING: low proton derivative',
     $  yprpmx, ' at bin',npmx,' of proton array'
        write(6,*) 'consider changing proton bg'
!         stop
        elseif(yprpmx.gt.1d16)then
        write(6,*) 'WARNING: high proton derivative',
     $  yprpmx, ' at bin',npmx,' of proton array'
        write(6,*) 'consider changing proton bg'
        endif 
        
        
        if(ypremx.lt.1d10) then
        write(6,*) 'WARNING: low electron derivative',
     $  ypremx, ' at bin',nemx,' of electron array'
        write(6,*) 'consider changing electron bg'
!         stop
        elseif(ypremx.gt.1d16)then
        write(6,*) 'WARNING: high electron derivative',
     $  ypremx, ' at bin',nemx,' of electron array'
        write(6,*) 'consider changing electron bg'        
        endif 
        
        if(yprgmx.lt.1d10) then
        write(6,*) 'WARNING: low photon derivative',
     $  ypremx, ' at bin',ngmx,' of photon array'
        write(6,*) 'consider changing photon bg'
!         stop
        elseif(yprgmx.gt.1d16)then
        write(6,*) 'WARNING: high photon derivative',
     $  yprgmx, ' at bin',ngmx,' of photon array'
        write(6,*) 'consider changing photon bg'        
        endif 
        
      endif    ! END of warnings for initial background choice    
       
	yprmax=abs(yprime(1))
	nmx=1

	do n=1,np+ne+ng+1+ne+np		!nt
	if (yprmax.lt.abs(yprime(n))) then
	yprmax=abs(yprime(n))
	nmx=n
	endif
        enddo

        sumptotl=0.
        
	do ik=1,10000

            if (ik*nout.eq.iprcl) then
	sumptotl=sumelpio+sumgrpio+sumntpio+sumnepio
     $  +summpio
	sumpioloss=sumptotl/sumprpio
C  Output on screen:
        write(6,1009) 
        write(6,*) 'step, time, bin of max deriv, max deriv'
	write (6,1400) iprcl,time,nmx,yprime(nmx)
	write(6,*) 'Thomson opt depth, photon compactness'
	write(6,1000) tauth, xden
	write(6,*) 'injection compactness: p,e,B'
	write (6,1000) exlumpr, exlumel,tb 
	write (6,*) 'total energy rates: p, e, g'
	write (6,1000) sumprot*ratmpe,sumelec,sumphot
	if(isyn.eq.1)then
	write (6,*) 'e-syn: loss, gain'
	write (6,1000) sumsyncel,sumsyncph
	endif
	if(iprsyn.eq.1)then
	write (6,*) 'p-syn: loss, gain'	
	write (6,1000) sumsyncpr*ratmpe, sumsyncprph
	endif
	if(icompt.eq.1)then
        write (6,*) 'e-ICS: Th loss, KN loss, tot loss, tot gain'
        write (6,1000) sumcsth,sumcskn,sumcsth+sumcskn, sumgcsth
        endif
        if(igg.eq.1)then
	write (6,*) 'gg->ee: loss, gain'	        
     	write (6,1000) sumggabs,sumgginj
     	endif 
     	if(ianni.eq.1)then 
        write (6,*) 'ee->g: loss, gain'
	write (6,1000) sumelann,sumanng
        endif 
! !  Output on file:     
            endif
        enddo 
	
1400  format(1x,i10,1x,1pe12.4,i5,1x,1pe12.4,1x,e12.4,1x,e12.4,1x,e12.4
     1    ,1x,e12.4)
1300  format(1x,1pe12.4,i5,1x,i5,1x,i5,2x,e12.4,2x,e12.4)
 
1000  format(1x,6(1pe12.4,1x))
1010  format('loss (%) via: e, g, ν, n, μ',1x,6(1pe12.4,1x))
1009  format(/)

	do 377 i=1,ntot
	   if( y(i).gt.1. .or. y(i).le.1.)then
              continue
           else
         print*,'overflow!!! i=',i,' y(i)=',y(i)
         call cpu_time(t_stop)
         write(*,*) 'Elapsed CPU time [s] = ', t_stop-t_start
         call itime(tarray)
         tt2=tarray(1)*3600.+tarray(2)*60.+tarray(3)
         write(*,*) 'Elapsed wallclock time [s] = ', tt2-tt1         
         stop
           end if
	   if( yprime(i).gt.1. .or. yprime(i).le.1.)then
              continue
           else
         print*,'overflow!!! i=',i,' yprime(i)=',yprime(i)
         call cpu_time(t_stop)
         write(*,*) 'Elapsed CPU time [s] = ', t_stop-t_start
         call itime(tarray)
         tt2=tarray(1)*3600.+tarray(2)*60.+tarray(3)
         write(*,*) 'Elapsed wallclock time [s] = ', tt2-tt1
         stop
           end if
377      continue

	tprev=time

      return
      end
************************************************************************

!!    PROTON INJECTION   
      subroutine prinj(gp,yp,np,deltap,extprinxp,sumprinj)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer ipexp
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION gp(np),yp(np),extprinxp(np)
     
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ap,iprext,ipexp
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr 

! sumprinj: energy injected into protons   
      sumprinj=0.
      if (iprext.eq.1) then
	 slpr=slprints
         q2=(gpextmx**(2.-slpr)-gpextmn**(2.-slpr))/(2.-slpr)
         qextpr=3.*exlumpr/q2
         do n=1,np 
        if(ipexp.eq.0) then
        if (gp(n).ge.gpextmn.and.gp(n).lt.gpextmx) then !if ipexp=0
        extprinxp(n)=qextpr/exp(yp(n))*gp(n)**(-slprints)
     $  *exp(-(gp(n)/gpextmx)**ap)**ipexp
         else
        extprinxp(n)=0.
        endif
        endif
        if(ipexp.eq.1) then
	if (gp(n).ge.gpextmn) then !if ipexp=1
        extprinxp(n)=qextpr/exp(yp(n))*gp(n)**(-slprints)
     $  *exp(-(gp(n)/gpextmx)**ap)**ipexp
        else
        extprinxp(n)=0.
        endif 
        endif
       sumprinj=sumprinj+deltap*gp(n)**2.*exp(yp(n))*extprinxp(n)
!        write(6,1000) gp(n), extprinxp(n), sumprinj, ap, ipexp
         enddo         
      endif
  
      if (iprextbr.eq.1) then
 	 slpr1=slprints1
	 slpr2=slprints2     
	 qp21=(gpextbr**(2.-slpr1)-gpextmn**(2.-slpr1))/(2.-slpr1)
         qp22=(gpextmx**(2.-slpr2)-gpextbr**(2.-slpr2))/(2.-slpr2)
	 fogp=gpextbr**(slpr1-slpr2)
	 qextpr=3.*exlumpr/(fogp*qp21+qp22)
	 do n=1,np
	 if (gp(n).ge.gpextmn.and.gp(n).le.gpextmx) then
	   if (gp(n).lt.gpextbr) then
	   extprinxp(n)=fogp*qextpr/exp(yp(n))*gp(n)**(-slpr1)
	   else
	 extprinxp(n)=qextpr/exp(yp(n))*gp(n)**(-slpr2)
     $   *exp(-(gp(n)/gpextmx)**ap2)**ipexp
	   end if  
         else
         extprinxp(n)=0.
         endif        
        sumprinj=sumprinj+deltap*gp(n)**2.*exp(yp(n))*extprinxp(n)
!         write(6,1000) gp(n), extprinxp(n), sumprinj
         enddo	
      endif 
      
      return 
1000  format(1x,6(1pe12.4,1x))      
      end
*************************************************************************

!!    ELECTRON INJECTION   
      subroutine elinj(ge,ye,ne,deltap,extelin,sumelinj)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer ieexp
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION ge(ne),ye(ne),extelin(ne)
      
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ae,ielext,ieexp  
      
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr

! sumelinj: energy injected into electrons    
      sumelinj=0.

      if (ielext.eq.1) then
         slel=slelints
         q2=(geextmx**(2.-slel)-geextmn**(2.-slel))/(2.-slel)
         qextel=3.*exlumel/q2
         do n=1,ne  
         if(ieexp.eq.0) then 
         if(ge(n).ge.geextmn.and.ge(n).le.geextmx)then !if ieexp=0
        extelin(n)=qextel/exp(ye(n))*ge(n)**(-slel)
     $  *exp(-(ge(n)/geextmx)**ae)**ieexp
         else
         extelin(n)=0.
         endif
         endif 
         if(ieexp.eq.1) then 
	if (ge(n).ge.geextmn) then !if ieexp=1
        extelin(n)=qextel/exp(ye(n))*ge(n)**(-slel)
     $  *exp(-(ge(n)/geextmx)**ae)**ieexp
         else
         extelin(n)=0.
         endif
         endif 
       sumelinj=sumelinj+deltap*ge(n)**2.*exp(ye(n))*extelin(n)
!         write(6,1000) log10(ge(n)), extelin(n)
         enddo         
      endif 
       
      if (ielextbr.eq.1) then
	 slel1=slelints1
	 slel2=slelints2    
         q21=(geextbr**(2.-slel1)-geextmn**(2.-slel1))/(2.-slel1)
         q22=(geextmx**(2.-slel2)-geextbr**(2.-slel2))/(2.-slel2)
         fog=geextbr**(slel1-slel2)
	 qextel=3.*exlumel/(fog*q21+q22)
	 do n=1,ne
	 if (ge(n).ge.geextmn.and.ge(n).le.geextmx) then
	   if (ge(n).lt.geextbr) then
	   extelin(n)=fog*qextel/exp(ye(n))*ge(n)**(-slel1)
	   else
	 extelin(n)=qextel/exp(ye(n))*ge(n)**(-slel2)
     $   *exp(-(ge(n)/geextmx)**ae2)**ieexp
	   end if  
         else
         extelin(n)=0.
         endif        
        sumelinj=sumelinj+deltap*ge(n)**2.*exp(ye(n))*extelin(n)
!        write(6,1000) ge(n), extelin(n), sumelinj
         enddo	
      endif 
      
      return 
1000  format(1x,6(1pe12.4,1x))   
      end       
    
*************************************************************************    

!!    ELECTRON SYNCHROTRON LOSSES
      subroutine synecool(ge,ye,syncel,sumsyncel)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION ge(ne),ye(ne),syncel(ne)
      
      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      
! sumsyncel: energy lost by electrons in synchrotron
	sumsyncel=0.

	do n=1,ne-1 
	syndf=(ge(n+1)**2.-1.)*exp(ye(n+1))-(ge(n)**2.-1.)*exp(ye(n))
	syncel(n)=4.*tb*syndf/deltap/ge(n)/exp(ye(n))/3.
        sumsyncel= sumsyncel+deltap*ge(n)**2.*syncel(n)*exp(ye(n))
        enddo
        return
      end         
*************************************************************************

!!    PROTON SYNCHROTRON LOSSES
      subroutine synpcool(gp,yp,syncpr,sumsyncpr)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION gp(np),yp(np),syncpr(np)

      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng

!  sumsyncpr: energy lost by protons in synchrotron
	sumsyncpr=0.

	do n=1,np-1
	syndf=(gp(n+1)**2.-1.)*exp(yp(n+1))-(gp(n)**2.-1.)*exp(yp(n))
	syncpr(n)=4.*tb*syndf/deltap/gp(n)/exp(yp(n))/3./ratmpe**3.
        sumsyncpr=sumsyncpr+deltap*gp(n)**2.*syncpr(n)*exp(yp(n))
        enddo
        return
      end
************************************************************************* 

!!    PHOTON ATTENUATION DUE TO SYNCHROTRON SELF ABSORPTION
      subroutine ssaphot(ge,ye,x,yg,sst, sumphssa)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,afine=7.3e-3,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION ge(ne),ye(ne),x(ng),yg(ng),sst(ng)
      
       common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 
     
!  uses delta-function approximation -- see equation (40) in Mastichiadis & Kirk (1995)
!  sumphssa: energy absorbed in synchrotron self absorption
        sumphssa=0.                
       do n=1,ng-1 
       if (ge(n).eq.0.) then
       sst(n)=0.
       else         
       elderiv=(exp(ye(n+1))/ge(n+1)**2.-exp(ye(n))/ge(n)**2.)/
     $ deltap
       if (elderiv.gt.0.) elderiv=-elderiv
       sst(n)=pi/6./afine/q**(.5)/x(n)**(.5)/ge(n)*elderiv
       endif
       if (x(n).gt.1.) sst(n)=0.
       sumphssa=sumphssa+xnorm*deltax*x(n)**2.*exp(yg(n))*sst(n)
       enddo  
        
        return
      end
*************************************************************************  

!!    ELECTRON INVERSE COMPTON LOSSES 
!     using 2 options for electron cooling
      subroutine compe_sim(x,ge,ye,yg,ikn,csth,cskn,sumcsth,
     $ sumcskn,taukn) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION ge(ne),ye(ne),x(ng),yg(ng), csth(ne), cskn(ne),
     $ gdens(ng),denkn(ng),p(npmax,npmax)
       
       common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 
       common/pyy/tauth 

       call densityold(x,ge,yg,ne,ng,deltax,gdens,denkn)
       
! sumcsth: electron energy losses due to ICS in Thomson regime
! sumcskn: electron energy losses due to ICS in KN regime
!  taukn: number of KN collisions which produce cool pairs

        sumcsth=0.
        sumcskn=0.
        taukn=0.
        if (ikn.eq.0) then
!  uses simplified scheme for KN losses -- see equation (45) in Mastichiadis & Kirk (1995)        
        do j=1, ne-1
	csdfth=(ge(j+1)**2.-1.)*exp(ye(j+1))-
     $	(ge(j)**2.-1.)*exp(ye(j))
	csth(j)=4.*tb*csdfth*q**(-5./3.)*gdens(j)/
     $          deltap/ge(j)/exp(ye(j))/3. 
        cskn(j)=xnorm*denkn(j)/ge(j)
        sumcsth=sumcsth+deltap*ge(j)**2.*exp(ye(j))*csth(j)
	sumcskn=sumcskn-deltap*ge(j)**2.*exp(ye(j))*cskn(j)
	taukn=taukn+deltap*cskn(j)*ge(j)*exp(ye(j))/tauth	
        enddo 
        endif 
                
        if(ikn.eq.1)then
!  uses full expression from Bloumenthal & Gould (1970)   
        call bgel(x,ge,yg,ye,ne,ng,deltap,deltax,p) 
        
        do j=1, ne-1
       
	sume1=0.
	sume2=0.

	do m=1,j
	if (m.eq.1.or.m.eq.j) then
	fcr=.5
	else
	fcr=1.
	endif
	fp1=fcr*p(j,m)*ge(m)*deltap
	sume1=sume1+fp1 
        enddo 

	do m=j+1,ne-1
	if (m.eq.1.or.m.eq.j) then
	fcr=.5
	else
	fcr=1.
	endif
	fp2=fcr*exp(ye(m))*ge(m)*p(m,j)*deltap
        sume2=sume2+fp2
        enddo

	csdf=(sume1-sume2/exp(ye(j)))
        cskn(j)=csdf 
        csdfth=(ge(j+1)**2.-1.)*exp(ye(j+1))-
     $	(ge(j)**2.-1.)*exp(ye(j))
	csth(j)=4.*tb*csdfth*q**(-5./3.)*gdens(j)/
     $          deltap/ge(j)/exp(ye(j))/3.  
        sumcsth=sumcsth+deltap*ge(j)**2.*exp(ye(j))*csth(j)
	sumcskn=sumcskn-deltap*ge(j)**2.*exp(ye(j))*cskn(j)
	taukn=taukn+deltap*cskn(j)*ge(j)*exp(ye(j))/tauth! to be checked
        enddo 
        endif 

        return
      end
************************************************************************* 
 
!!    PHOTON LOSSES DUE TO PHOTON-PHOTON PAIR PRODUCTION
      subroutine ggphot(x,yg,ggabsr,sumggabs) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION x(ng), yg(ng), ggabsr(ng), yggal(ng)

       common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 
       
        call ggabsal(deltax,x,yg,ng,yggal)
        
! sumggabs: energy lost due to photon-photon absortion
        sumggabs=0.
	do n=1,ng
        ggabsr(n)=tb*q**(-5./3.)*yggal(n)
	sumggabs=sumggabs+deltax*x(n)**2.*ggabsr(n)*exp(yg(n))*xnorm
        enddo 
        return
      end
*************************************************************************    

!!    PHOTON INJECTION FROM ELECTRON SYNCHROTRON
      subroutine synephot(ge,ye,x,yg,syncph,sumsyncph)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION ge(ne),ye(ne),x(ng),yg(ng),syncph(ng),sumsy(ng)
      
      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
	
! sumsyncph: energy injected in photons by electron synchrotron losses
	sumsyncph=0.

        do n=1,ng
	sumsy(n)=0.
! for each photon energy we form the ratio x/xcrit
            do m=1,ne
! 	xsycr=1.5*ge(m)**2.*q
	xsycr=1.5*(ge(m)**2.-1.)*q !19/11/20
	ysycr=x(n)/xsycr
            if (ysycr.gt.30.) then
	fsyval=0.
	  else
	fsyval=fsynch(ysycr)/4./pi
            end if
	sumsy(n)=sumsy(n)+deltap*ge(m)*exp(ye(m))*
     $  (4.*sqrt(3.)*tb/q/x(n))*fsyval/xnorm/exp(yg(n))
            enddo 
        enddo

	do n=1,ng
	syncph(n)=sumsy(n)
	sumsyncph = sumsyncph +
     $    deltax*xnorm*x(n)**2.*exp(yg(n))*syncph(n)	
        enddo
        return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************

!!    PHOTON INJECTION FROM PROTON SYNCHROTRON
      subroutine synpphot(gp,yp,x,yg,syncprph,sumsyncprph)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION gp(np),yp(np),x(ng),yg(ng),syncprph(ng),sumpsy(ng)
      
      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
 
!  sumsyncprph: energy injected in photons by proton synchrotron losses
	sumsyncprph=0.

	do n=1,ng
	sumpsy(n)=0.
! for each photon energy we form the ratio x/xcrit
            do m=1,np
	xsycr=1.5*gp(m)**2.*q/ratmpe
	ysycr=x(n)/xsycr
	if (ysycr.gt.30.) then
	fsyval=0.
	  else
	fsyval=fsynch(ysycr)/4./pi
	end if
	sumpsy(n)=sumpsy(n)+deltap*gp(m)*exp(yp(m))*
     $      (4.*sqrt(3.)*tb/q/x(n))*fsyval/xnorm/exp(yg(n))
	   enddo
        enddo 

	do n=1,ng
	syncprph(n)=sumpsy(n)/ratmpe
	sumsyncprph = sumsyncprph +
     $    deltax*xnorm*x(n)**2.*exp(yg(n))*syncprph(n)
        enddo 
        
        return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************

!!    PHOTON INJECTION FROM ELECTRON INVERSE COMPTON SCATTERING
      subroutine compph_bg(x,ge,yg,ye,gcsth, sumgcsth)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION  ge(ne),ye(ne),x(ng),yg(ng),gcsth(ng)
      
      common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
      
      call bg(x,ge,yg,ye,ne,ng,deltap,deltax,gcsth) 
! sumgcsth: energy injected in photons from electron inverse Compton scattering
       sumgcsth=0.
       do j=1,ng
       sumgcsth=sumgcsth+deltax*x(j)**2.*gcsth(j)*exp(yg(j))*xnorm
       enddo 

       return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************
 
!!    PAIR INJECTION DUE TO PHOTON-PHOTON PAIR PRODUCTION
      subroutine ggelec(x,yg,ge,ye,gginje,sumggabs,sumgginj) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION x(ng),yg(ng),ge(ne),ye(ne),gginje(ne),yggcral(ng)

       common/param/tb,q,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 

        call ggcreal(deltax,x,ge,yg,ne,ng,yggcral)
        
! sumgginj: energy injection into pairs due to photon-photon pair production
        sumgginj=0.
        tempgginj=0.
	do n=1,ne
        gginje(n)=4.*tb**2.*q**(-10./3.)*yggcral(n)/exp(ye(n))
	tempgginj=tempgginj+deltap*ge(n)**2.*gginje(n)*exp(ye(n))
        enddo
                
! re-normalization based on photon energy losses 
        do n=1,ne
        gginje(n)=gginje(n)*sumggabs/tempgginj
        sumgginj=sumgginj+deltap*ge(n)**2.*gginje(n)*exp(ye(n))
        enddo 

        return        
       end 
*************************************************************************

      subroutine densityold(x,ge,yg,ne,ng,deltax,gdens,denkn)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),gdens(ng),denkn(ng)

      do 1 n=1,ne,2
      xtrns=3./4./ge(n)
        
        do m=1,ng
	  if (x(m).gt.xtrns) then
	  mtr=m-1
	  goto 3
	  endif
        enddo

	
3	sumth=.0
	if (mtr.le.1) goto 8

	do m=1,mtr
	gx=x(m)*ge(n)
! 	fkn=1.-gx
	fkn=1.	
	  if (m.eq.1.or.m.eq.mtr) then
	     fac=.5
 	  else
	     fac=1.
	  endif
	sumth=sumth+fac*fkn*deltax*x(m)**2.*exp(yg(m))
        enddo
   
8	sumkn=0.
	do m=mtr+1,ng
	  if (m.eq.mtr+1.or.m.eq.ng) then
	     fac=.5
 	  else
	     fac=1.
	  endif
	ylog=log(4.*x(m)*ge(n)-11./6.)
	sumkn=sumkn+fac*deltax*exp(yg(m))*ylog
        enddo 

        gdens(n)=sumth
        if (gdens(n).eq.0.) gdens(n)=1.e-40
	denkn(n)=sumkn
!   	write (6,1000) n,ge(n),gdens(n)
1       continue

	do n=2,ne-1,2

	if (gdens(n-1).gt.0.) then 
	gdlm=log(gdens(n-1))
	gdlp=log(gdens(n+1))
	gdl=gdlm+.5*(gdlp-gdlm)
        gdens(n)=exp(gdl)
	endif

	dknm=log(denkn(n-1))
	dknp=log(denkn(n+1))
	dkn=dknm+.5*(dknp-dknm)
	denkn(n)=exp(dkn)
        enddo

1000    format (2x,i10,2x,4(1pd12.4,2x))
	end
*************************************************************************

       subroutine bg(x,ge,yg,ye,ne,ng,deltap,deltax,gcsth)
       implicit real*8(a-h,o-z)
       dimension x(ng),ge(ne),yg(ng),ye(ne),gcsth(ng)

       do l=1,ng
        xsc=x(l)

	sumg=0.
	do 15 m=1,ne
	if (ge(m).lt.xsc) then
	goto 15
	endif

	sumx=0.
	do 20 n=1,ng
	xsoft=x(n)	

	if (xsoft.ge.xsc) then
	goto 20
	endif

	if (xsoft.gt.ge(m)) then
	goto 20
	endif

	e0mn=xsc/4./(ge(m)-xsc)/ge(m)
	if (xsoft.le.e0mn) goto 20
	
	Esc=xsc/ge(m)
	Gebg=4.*ge(m)*xsoft

        qs = Esc/(Gebg*(1.-Esc))
	
        sa=2.*qs*log(qs)+(1.+2.*qs)*(1.-qs)
        sb=(gebg*qs)**2.*(1.-qs)/(2.*(1+Gebg*qs))
	Fg=sa+sb
	
	if (Fg.lt.0.) then
	goto 20
	endif

	fux=deltax*exp(yg(n))*Fg
	sumx=sumx+fux
	
	if (sumx.gt.0.) then
	if (fux/sumx.lt.1.e-4) goto 21
	endif
20 	continue

21	fug=deltap*sumx*exp(ye(m))/ge(m)
	if (sumg.gt.0.) then
	if (fug/sumg.lt.1.e-4) goto 16
	endif
14	sumg=sumg+fug

15	continue
	
16	continue
	gcsth(l)=.75*sumg/exp(yg(l))
        enddo 
        
1000    format (1x,5(1pd12.4,1x))
	return
	end 
*************************************************************************

      subroutine bgel(x,ge,yg,ye,ne,ng,deltap,deltax,p)
! we construct the probability function P(ga,gb)--see equation (5.17) in Bloumenthal & Gould (1970)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),ye(ne),p(npmax,npmax)
      dimension xa(2*ng),yga(2*ng)
      common/tfl/xnorm

! we divide the photons into 2*ng bins
	do 21 m=1,2*ng-1
21	xa(m)=x(1)*exp(deltap*(m-1)) 

	do 22 m=1,2*ng-1,2 
	n=(m+1)/2
22      yga(m)=yg(n)

	do 23 m=2,2*ng-2,2
	n=m/2 
23      yga(m)=.5*(yg(n)+yg(n+1))


	do 10 k=1,ne
	ga=ge(k)

	do 11 l=1,k
	gb=ge(l)

	if (k.eq.l) p(k,l)=0.
 
	xsc=ga-gb
	xmn=4./gb/8. !May 2017 -- original xmn=4./gb
	sumx=0.	 				
	
	do 12 m=1,2*ng-1
	if (xa(m).lt.xmn) then
	kmin=m
	goto 12
	endif

	if (ge(k).le.1.*xa(m)) goto110
	
	if (xa(m).ge.xsc) goto110

	e0mn=xa(m)/4./(ge(k)-xsc)/ge(k)
	e0mn=xsc/4./(ge(k)-xsc)/ge(k)
	if (xa(m).lt.e0mn) goto12
	
	Esc=xsc/ge(k)
	Gebg=4.*ge(k)*xa(m)

        qs = Esc/(Gebg*(1.-Esc))
	
        sa=2.*qs*log(qs)+(1.+2.*qs)*(1.-qs)
        sb=(gebg*qs)**2.*(1.-qs)/(2.*(1+Gebg*qs))
	Fg=sa+sb
	
	if (Fg.lt.0.) goto 110

	fux=deltap*exp(yga(m))*xnorm*Fg
	sumx=sumx+fux
!       stop integration when the relative change in sumx is less than 1e-4      
!!      CAUTION: this choice makes the code faster, but it may give wrong results when deep in KN. 
!!      In that case, comment out the following if statement
	if (sumx.gt.0.) then   
	if (fux/sumx.lt.1.e-4) goto 110
	endif
12 	continue
110	p(k,l)=.75*sumx/ge(k)**2.
11	continue

10	continue
	
1000    format (1x,5(1pd12.4,1x))
	return
	end 
*************************************************************************

        subroutine ggabsal(deltax,x,yg,ng,ygi)
        implicit real*8(a-h,o-z)
        dimension x(ng),yg(ng),ygi(ng)
!       uses the approximate expression from Coppi & Blandford (1990)
!       see also equation (55) in Mastichiadis & Kirk (1990)
	xmin=x(1)
	xmax=x(ng)

	do n=1,ng
	xinv=1./x(n)
	sum=0.
            do m=1,ng
	if (x(m).gt.xinv) then
	w=x(m)*x(n)
	r=0.652*(w**2.-1.)*log(w)/w**3.
	sum=sum+deltax*r*exp(yg(m))*x(m)
	end if
            enddo
	ygi(n)=sum
        enddo 
        return
	end
*************************************************************************

      subroutine ggcreal(deltax,x,ge,yg,ne,ng,yggcre)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),yga(500),yggcre(ng)

	do 1 n=1,ne
	sum=0.
	if (ge(n)**2.lt.2.) then
        yggcre(n)=0.
	goto 1
	endif
	
	xa=ge(n)+sqrt(ge(n)**2-2.)
	xa=2.*ge(n)
	xb=1./xa

	if (xa.gt.x(ng).or.xa.lt.x(1)) then
	  yga(n)=-99.
	  goto 12
	end if

	if (xb.gt.x(ng).or.xb.lt.x(1)) then
	  yga(n)=-99.
	  goto 12
	end if

	do ka=1,ng
	 if (x(ka).gt.xa) then
	   ma=ka-1
	   goto 30
	 endif
        enddo

30        yga(n)=yg(ma+1)+(yg(ma)-yg(ma+1))*log(xa/x(ma+1))
     $    /log(x(ma)/x(ma+1))

	do m=1,ng
	if (x(m).gt.xb) then
	w=x(m)*xa
	r=0.652*(w**2.-1.)*log(w)/w**3.
	sum=sum+deltax*r*exp(yg(m))*x(m)
	end if
        enddo

12      yggcre(n)=exp(yga(n))*sum
!  	write (6,1000) ge(n),yggcre(n)

1       continue
1000    format (2x,4(1pd12.4,2x))
	end
*************************************************************************
  
      
************************************************************************
      FUNCTION FSYNCH(T)
! 
!  THE FUNCTION F(X)=X * INTEGRAL(X TO INFINITY) K_5/3(Z) DZ
!  WHICH IS USED IN THE CALCULATION OF SYNCHROTRON EMISSIVITY.
!  (SEE  'RADIO ASTROPHYSICS' BY A.G.PACHOLCZYK, 1970, (FREEMAN,
!  SAN FRANCISCO).
!                                 R.J.PROTHEROE, ADELAIDE 20 OCT 1988.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X1(35),Y1(35),C1(34,3)
      DIMENSION X2(126),Y2(126),C2(125,3),C21(125),C22(125),C23(125)
      EQUIVALENCE (C2(1,1),C21(1)),(C2(1,2),C22(1)),(C2(1,3),C23(1))
      F1(X,A1)=A1*(X/2.D0)**(1.D0/3.D0)
      F2(X,A2)=A2*EXP(0.5D0*DLOG(X)-X)
!
      DATA A1/2.7082359D0/,A2/1.25331414/
      DATA X1/
     *  1.00D-4, 2.00D-4, 5.00D-4, 1.00D-3, 2.00D-3, 5.00D-3, 1.00D-2,
     *  2.00D-2, 3.00D-2, 4.00D-2, 5.00D-2, 6.00D-2, 7.00D-2, 8.00D-2,
     *  9.00D-2, 1.00D-1, 1.10D-1, 1.20D-1, 1.30D-1, 1.40D-1, 1.50D-1,
     *  1.60D-1, 1.70D-1, 1.80D-1, 1.90D-1, 2.00D-1, 2.10D-1, 2.20D-1,
     *  2.30D-1, 2.40D-1, 2.50D-1, 2.60D-1, 2.70D-1, 2.80D-1, 2.90D-1/
      DATA X2/
     *  2.90D-1, 3.00D-1, 3.10D-1, 3.20D-1, 3.30D-1, 3.40D-1, 3.50D-1,
     *  3.60D-1, 3.70D-1, 3.80D-1, 3.90D-1, 4.00D-1, 4.10D-1, 4.20D-1,
     *  4.30D-1, 4.40D-1, 4.50D-1, 4.60D-1, 4.70D-1, 4.80D-1, 4.90D-1,
     *  5.00D-1, 5.20D-1, 5.40D-1, 5.60D-1, 5.80D-1, 6.00D-1, 6.20D-1,
     *  6.40D-1, 6.60D-1, 6.80D-1, 7.00D-1, 7.20D-1, 7.40D-1, 7.60D-1,
     *  7.80D-1, 8.00D-1, 8.20D-1, 8.40D-1, 8.60D-1, 8.80D-1, 9.00D-1,
     *  9.20D-1, 9.40D-1, 9.60D-1, 9.80D-1, 1.00D+0, 1.05D+0, 1.10D+0,
     *  1.15D+0, 1.20D+0, 1.25D+0, 1.30D+0, 1.35D+0, 1.40D+0, 1.45D+0,
     *  1.50D+0, 1.55D+0, 1.60D+0, 1.65D+0, 1.70D+0, 1.75D+0, 1.80D+0,
     *  1.85D+0, 1.90D+0, 1.95D+0, 2.00D+0, 2.10D+0, 2.20D+0, 2.30D+0,
     *  2.40D+0, 2.50D+0, 2.60D+0, 2.70D+0, 2.80D+0, 2.90D+0, 3.00D+0,
     *  3.10D+0, 3.20D+0, 3.30D+0, 3.40D+0, 3.50D+0, 3.60D+0, 3.70D+0,
     *  3.80D+0, 3.90D+0, 4.00D+0, 4.10D+0, 4.20D+0, 4.30D+0, 4.40D+0,
     *  4.50D+0, 4.60D+0, 4.70D+0, 4.80D+0, 4.90D+0, 5.00D+0, 5.25D+0,
     *  5.50D+0, 5.75D+0, 6.00D+0, 6.25D+0, 6.50D+0, 6.75D+0, 7.00D+0,
     *  7.25D+0, 7.50D+0, 7.75D+0, 8.00D+0, 8.25D+0, 8.50D+0, 8.75D+0,
     *  9.00D+0, 9.25D+0, 9.50D+0, 9.75D+0, 1.00D+1, 1.20D+1, 1.40D+1,
     *  1.60D+1, 1.80D+1, 2.00D+1, 2.50D+1, 3.00D+1, 4.00D+1, 5.00D+1/
      DATA Y1/
     *  1.00D+0, 9.94D-1, 9.96D-1, 9.91D-1, 9.97D-1, 9.74D-1, 9.61D-1,
     *  9.37D-1, 9.18D-1, 9.02D-1, 8.86D-1, 8.71D-1, 8.58D-1, 8.44D-1,
     *  8.32D-1, 8.20D-1, 8.09D-1, 7.97D-1, 7.86D-1, 7.76D-1, 7.65D-1,
     *  7.56D-1, 7.47D-1, 7.37D-1, 7.28D-1, 7.19D-1, 7.11D-1, 7.02D-1,
     *  6.93D-1, 6.85D-1, 6.77D-1, 6.68D-1, 6.60D-1, 6.53D-1, 6.45D-1/
      DATA Y2/
     *  1.82D+0, 1.80D+0, 1.79D+0, 1.78D+0, 1.77D+0, 1.76D+0, 1.75D+0,
     *  1.73D+0, 1.72D+0, 1.71D+0, 1.70D+0, 1.70D+0, 1.69D+0, 1.68D+0,
     *  1.67D+0, 1.66D+0, 1.66D+0, 1.65D+0, 1.64D+0, 1.63D+0, 1.63D+0,
     *  1.62D+0, 1.61D+0, 1.59D+0, 1.58D+0, 1.57D+0, 1.56D+0, 1.55D+0,
     *  1.54D+0, 1.53D+0, 1.52D+0, 1.51D+0, 1.50D+0, 1.50D+0, 1.49D+0,
     *  1.48D+0, 1.47D+0, 1.46D+0, 1.46D+0, 1.45D+0, 1.44D+0, 1.44D+0,
     *  1.43D+0, 1.43D+0, 1.42D+0, 1.42D+0, 1.42D+0, 1.40D+0, 1.39D+0,
     *  1.38D+0, 1.37D+0, 1.36D+0, 1.35D+0, 1.34D+0, 1.33D+0, 1.32D+0,
     *  1.30D+0, 1.30D+0, 1.29D+0, 1.29D+0, 1.29D+0, 1.28D+0, 1.27D+0,
     *  1.27D+0, 1.27D+0, 1.27D+0, 1.25D+0, 1.23D+0, 1.23D+0, 1.23D+0,
     *  1.23D+0, 1.23D+0, 1.23D+0, 1.21D+0, 1.20D+0, 1.20D+0, 1.20D+0,
     *  1.20D+0, 1.18D+0, 1.19D+0, 1.20D+0, 1.19D+0, 1.19D+0, 1.19D+0,
     *  1.19D+0, 1.19D+0, 1.18D+0, 1.19D+0, 1.19D+0, 1.19D+0, 1.17D+0,
     *  1.15D+0, 1.17D+0, 1.14D+0, 1.16D+0, 1.13D+0, 1.13D+0, 1.13D+0,
     *  1.12D+0, 1.12D+0, 1.11D+0, 1.10D+0, 1.09D+0, 1.09D+0, 1.09D+0,
     *  1.08D+0, 1.09D+0, 1.08D+0, 1.08D+0, 1.10D+0, 1.08D+0, 1.07D+0,
     *  1.07D+0, 1.07D+0, 1.07D+0, 1.07D+0, 1.07D+0, 1.06D+0, 1.05D+0,
     *  1.04D+0, 1.04D+0, 1.04D+0, 1.03D+0, 1.02D+0, 1.02D+0, 1.00D+0/
      DATA C1/
     * -8.88D-16,-6.35D+1, 2.34D+1,-1.78D+1, 1.22D+1,-1.26D+1, 9.44D-1,
     * -3.09D+0,-1.50D+0,-1.58D+0,-1.56D+0,-1.43D+0,-1.31D+0,-1.35D+0,
     * -1.22D+0,-1.11D+0,-1.15D+0,-1.16D+0,-1.03D+0,-1.05D+0,-1.01D+0,
     * -9.17D-1,-9.23D-1,-9.09D-1,-9.24D-1,-8.81D-1,-8.41D-1,-8.83D-1,
     * -8.58D-1,-8.08D-1,-8.34D-1,-8.27D-1,-7.84D-1,-7.58D-1,-1.05D+6,
     *  4.13D+5,-1.24D+5, 4.16D+4,-1.16D+4, 3.30D+3,-5.82D+2, 1.78D+2,
     * -1.87D+1, 1.06D+1,-8.35D+0, 2.19D+1,-1.07D+1, 7.46D+0, 5.25D+0,
     *  5.26D+0,-8.76D+0, 7.92D+0, 4.74D+0,-7.09D+0, 1.12D+1,-1.64D+0,
     *  1.03D+0, 3.88D-1,-1.96D+0, 6.25D+0,-2.23D+0,-1.91D+0, 4.37D+0,
     *  5.97D-1,-3.21D+0, 3.95D+0, 3.28D-1, 2.29D+0, 4.87D+9,-5.97D+8,
     *  1.10D+8,-1.77D+7, 1.65D+6,-2.59D+5, 2.53D+4,-6.56D+3, 9.76D+2,
     * -6.31D+2, 1.01D+3,-1.09D+3, 6.05D+2,-7.34D+1, 7.24D-2,-4.67D+2,
     *  5.56D+2,-1.06D+2,-3.94D+2, 6.10D+2,-4.28D+2, 8.90D+1,-2.14D+1,
     * -7.81D+1, 2.74D+2,-2.83D+2, 1.05D+1, 2.10D+2,-1.26D+2,-1.27D+2,
     *  2.39D+2,-1.21D+2, 6.55D+1,-9.82D+1/
      DATA C21/
     * -1.32D+0,-1.30D+0,-1.29D+0,-1.20D+0,-1.08D+0,-1.12D+0,-1.15D+0,
     * -1.06D+0,-9.53D-1,-1.05D+0,-9.44D-1,-8.96D-1,-8.40D-1,-7.85D-1,
     * -7.50D-1,-6.57D-1,-7.86D-1,-6.62D-1,-7.83D-1,-7.48D-1,-5.46D-1,
     * -6.18D-1,-7.86D-1,-6.19D-1,-5.08D-1,-5.38D-1,-5.33D-1,-5.57D-1,
     * -5.27D-1,-4.25D-1,-3.95D-1,-4.52D-1,-3.89D-1,-4.52D-1,-4.37D-1,
     * -3.27D-1,-3.92D-1,-4.53D-1,-3.69D-1,-2.79D-1,-3.64D-1,-3.74D-1,
     * -2.08D-1,-2.01D-1,-4.12D-2,-4.82D-3,-2.72D-1,-2.72D-1,-2.32D-1,
     * -2.85D-1,-1.50D-1,-2.74D-1,-1.50D-1,-1.47D-1,-2.50D-1,-2.82D-1,
     * -2.30D-1,-4.62D-2,-1.09D-1, 9.37D-3,-1.56D-1,-1.79D-1,-9.21D-2,
     *  1.84D-2,-2.42D-3,-2.14D-1,-2.50D-1,-1.20D-1, 7.88D-2,-2.46D-2,
     * -1.51D-2,-2.43D-2,-7.13D-2,-1.61D-1,-6.90D-2, 4.10D-2,-1.40D-3,
     * -1.43D-1,-5.99D-2, 1.24D-1, 1.27D-2,-2.62D-2, 1.83D-2,-1.29D-2,
     * -7.50D-3,-1.06D-1,-5.70D-3, 1.06D-1,-1.84D-2,-5.20D-2,-2.98D-1,
     *  3.12D-2,-7.91D-2,-3.64D-2, 3.25D-2,-2.25D-1, 8.50D-2,-4.70D-2,
     * -1.03D-2,-2.57D-2,-4.80D-2,-3.25D-2,-7.97D-3,-1.54D-2,-1.28D-2,
     * -8.27D-3, 1.09D-2,-3.91D-2, 4.01D-2, 5.86D-3,-6.76D-2, 1.17D-2,
     * -1.24D-2,-1.12D-2, 1.96D-3,-2.00D-3,-1.31D-2,-2.62D-3,-3.72D-3,
     * -2.18D-3,-1.30D-3,-5.91D-4,-1.91D-3,-4.06D-4,-1.72D-3/
      DATA C22/
     *  4.06D+0,-2.76D+0, 4.47D+0, 3.79D+0, 8.80D+0,-1.26D+1, 8.78D+0,
     *  8.23D-1, 9.69D+0,-1.92D+1, 2.97D+1,-2.50D+1, 3.06D+1,-2.51D+1,
     *  2.87D+1,-1.95D+1, 6.59D+0, 5.84D+0,-1.79D+1, 2.15D+1,-1.36D+0,
     * -5.82D+0,-2.59D+0, 1.10D+1,-5.42D+0, 3.89D+0,-3.64D+0, 2.45D+0,
     * -9.30D-1, 6.01D+0,-4.49D+0, 1.61D+0, 1.57D+0,-4.71D+0, 5.46D+0,
     *  3.99D-2,-3.30D+0, 2.36D-1, 3.95D+0, 5.60D-1,-4.79D+0, 4.30D+0,
     *  4.01D+0,-3.66D+0, 1.16D+1,-9.82D+0,-3.55D+0, 3.56D+0,-2.77D+0,
     *  1.69D+0, 1.01D+0,-3.49D+0, 5.96D+0,-5.91D+0, 3.86D+0,-4.50D+0,
     *  5.55D+0,-1.87D+0, 6.19D-1, 1.75D+0,-5.05D+0, 4.60D+0,-2.86D+0,
     *  5.07D+0,-5.49D+0, 1.27D+0,-2.00D+0, 3.31D+0,-1.32D+0, 2.86D-1,
     * -1.90D-1, 9.75D-2,-5.67D-1,-3.30D-1, 1.25D+0,-1.49D-1,-2.75D-1,
     * -1.14D+0, 1.97D+0,-1.24D-1,-9.91D-1, 6.02D-1,-1.57D-1,-1.56D-1,
     *  2.10D-1,-1.19D+0, 2.20D+0,-1.08D+0,-1.63D-1,-1.73D-1,-2.29D+0,
     *  5.58D+0,-6.68D+0, 7.11D+0,-6.42D+0, 3.84D+0,-7.40D-1, 2.12D-1,
     * -6.59D-2, 4.31D-3,-9.35D-2, 1.56D-1,-5.74D-2, 2.76D-2,-1.70D-2,
     *  3.51D-2, 4.17D-2,-2.42D-1, 5.58D-1,-6.95D-1, 4.01D-1,-8.43D-2,
     * -1.21D-2, 1.71D-2, 3.54D-2,-5.12D-2, 6.95D-3,-1.72D-3, 1.17D-3,
     * -3.97D-4, 8.37D-4,-4.84D-4, 2.20D-4, 8.02D-5,-2.12D-4/
      DATA C23/
     * -2.27D+2, 2.41D+2,-2.29D+1, 1.67D+2,-7.14D+2, 7.13D+2,-2.65D+2,
     *  2.96D+2,-9.64D+2, 1.63D+3,-1.82D+3, 1.85D+3,-1.86D+3, 1.80D+3,
     * -1.61D+3, 8.69D+2,-2.48D+1,-7.93D+2, 1.31D+3,-7.61D+2,-1.49D+2,
     *  5.39D+1, 2.26D+2,-2.73D+2, 1.55D+2,-1.26D+2, 1.02D+2,-5.64D+1,
     *  1.16D+2,-1.75D+2, 1.02D+2,-7.04D-1,-1.05D+2, 1.69D+2,-9.03D+1,
     * -5.57D+1, 5.90D+1, 6.20D+1,-5.66D+1,-8.92D+1, 1.51D+2,-4.88D+0,
     * -1.28D+2, 2.55D+2,-3.58D+2, 1.04D+2, 4.74D+1,-4.22D+1, 2.97D+1,
     * -4.55D+0,-3.00D+1, 6.30D+1,-7.91D+1, 6.51D+1,-5.57D+1, 6.70D+1,
     * -4.95D+1, 1.66D+1, 7.52D+0,-4.53D+1, 6.43D+1,-4.97D+1, 5.29D+1,
     * -7.04D+1, 4.51D+1,-2.18D+1, 1.77D+1,-1.54D+1, 5.35D+0,-1.59D+0,
     *  9.59D-1,-2.21D+0, 7.89D-1, 5.27D+0,-4.66D+0,-4.20D-1,-2.88D+0,
     *  1.03D+1,-6.97D+0,-2.89D+0, 5.31D+0,-2.53D+0, 1.66D-3, 1.22D+0,
     * -4.68D+0, 1.13D+1,-1.09D+1, 3.05D+0,-3.17D-2,-7.04D+0, 2.62D+1,
     * -4.08D+1, 4.60D+1,-4.51D+1, 3.42D+1,-1.53D+1, 1.27D+0,-3.71D-1,
     *  9.37D-2,-1.30D-1, 3.32D-1,-2.84D-1, 1.13D-1,-5.94D-2, 6.95D-2,
     *  8.79D-3,-3.78D-1, 1.07D+0,-1.67D+0, 1.46D+0,-6.47D-1, 9.63D-2,
     *  3.90D-2, 2.43D-2,-1.16D-1, 7.76D-2,-1.45D-3, 4.82D-4,-2.61D-4,
     *  2.06D-4,-2.20D-4, 4.69D-5,-9.34D-6,-9.74D-6, 1.99D-5/
C
      IF(T.LE.0.D0)GOTO 103
      IF(T.LT.X1(1))GOTO 101
      IF(T.GE.X2(126))GOTO 102
      IF(T.GE.X1(35))GOTO 20
C
      I1=1
      IF(T.GE.X1(9))I1=9
      IF(T.GE.X1(18))I1=18
      IF(T.GE.X1(27))I1=27
      DO 10 K=I1,(I1+8)
      I=K
      IF(X1(I+1).GT.T)GOTO 11
   10 CONTINUE
   11 D=T-X1(I)
      S=((C1(I,3)*D+C1(I,2))*D+C1(I,1))*D+Y1(I)
      FSYNCH=F1(T,A1)*S
      RETURN
C
   20 I1=1
      IF(T.GE.X2(10))I1=10
      IF(T.GE.X2(20))I1=20
      IF(T.GE.X2(30))I1=30
      IF(T.GE.X2(40))I1=40
      IF(T.GE.X2(50))I1=50
      IF(T.GE.X2(60))I1=60
      IF(T.GE.X2(70))I1=70
      IF(T.GE.X2(80))I1=80
      IF(T.GE.X2(90))I1=90
      IF(T.GE.X2(100))I1=100
      IF(T.GE.X2(110))I1=110
      IF(T.GE.X2(120))I1=120
      DO 22 K=I1,(I1+9)
      I=K
      IF(X2(I+1).GT.T)GOTO 21
   22 CONTINUE
   21 D=T-X2(I)
      S=((C2(I,3)*D+C2(I,2))*D+C2(I,1))*D+Y2(I)
      FSYNCH=F2(T,A2)*S
      RETURN
  101 FSYNCH=F1(T,A1)
      RETURN
  102 FSYNCH=F2(T,A2)
      RETURN
  103 FSYNCH=0.D0
      RETURN
      END



*************************************************************************
      