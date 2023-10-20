      subroutine mc_shms_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,
     > wcs_flag,
     > decay_flag,dflag,resmult,spec,ok_hut,zinit,pathlen,
     > spectr,xs_pad1x_x,xs_pad1x_y,ys_pad1y_x,ys_pad1y_y,
     > xs_pad2x_x,xs_pad2x_y,ys_pad2y_x,ys_pad2y_y)
C----------------------------------------------------------------------
C
C Monte-Carlo of HMS detector hut.
C----------------------------------------------------------------------

	implicit 	none
	include 'struct_shms.inc'
	include '../spectrometers.inc'
	include 'hut.inc'

C Math constants
	integer hit_s1x,hit_s1y,hit_s2x,hit_s2y
	real*8 pi,d_r,r_d,root

	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2

	logical*4 cer_flag
	logical*4 vac_flag
        parameter (cer_flag = .false.) ! TRUE means 1st Cerenkov (Ar/Ne) is in front of chambers
        parameter (vac_flag = .true.) ! FALSE means helium bag replaces 1st Cerenkov (Ar/Ne) 

c	common /hutflag/ cer_flag,vac_flag
	common /paddles_x/pad_1x_lo_num,pad_1x_hi_num,pad_2x_lo_num
     >  ,pad_2x_hi_num
        common /paddles_y/pad_1y_lo_num,pad_1y_hi_num,pad_2y_lo_num
     >  ,pad_2y_hi_num

C all parameters, later to take from .parm files
C The arguments

	real*8	p,m2			!momentum and mass of particle
	real*8	x_fp,y_fp,dx_fp,dy_fp	!Focal plane values to return
	logical	ms_flag, wcs_flag	!particle, m_scat, wc_smear
	logical	ok_hut			!true if particle makes it
	logical decay_flag,dflag
	real*8	xcal,ycal		!Position of track at calorimeter.
	real*8 pathlen
	real*8 resmult				!DC resolution factor
	real*8 zinit
        real*4 spec(58)
	real*8 xs_pad1x_x,xs_pad1x_y
        real*8 xs_pad2x_x,xs_pad2x_y
        real*8 ys_pad1y_x,ys_pad1y_y
	real*8 ys_pad2y_x,ys_pad2y_y

c external function
	real*8 gauss1
C Local declarations.

	integer*4	spectr	!spectrometer number (for tune-dependent stuff)
	integer*4	i,iplane,jchamber,npl_off
	integer*4	chan	/1/

	real*8	radw,drift
c need real*4 for cernlib routine lfut
	real*4	xdc(12),ydc(12),zdc(12)		!positions at d.c. planes
	real*4	x0,y0,dx,dy			!track fitting temporaries
	real*4	badf				!temporaries

	real*8  tmpran1,tmpran2
	real*8   nsig_max
	real*8 hdc_del_plane
	parameter (nsig_max=99.0d0)
c mkj
	logical use_det_cut
	parameter (use_det_cut=.true.)
!	parameter (use_det_cut=.false.)
! cadd cut on the hodoscope paddle                                                                                                                                                                                                                                            \                                                                                                                                                                                                                                   
        real*8 pscin_1x_center,pscin_1x_size,pscin_1x_spacing
        real*8 pscin_2x_center,pscin_2x_size,pscin_2x_spacing
        real*8 pscin_1y_center,pscin_1y_size,pscin_1y_spacing
        real*8 pscin_2y_center,pscin_2y_size,pscin_2y_spacing
        real*8 pscin_2y_size_hi
        real*8 xs_lo,xs_hi,ys_lo,ys_hi
        integer pad_1x_lo_num,pad_1x_hi_num,pad_2x_lo_num,pad_2x_hi_num
        integer pad_1y_lo_num,pad_1y_hi_num,pad_2y_lo_num,pad_2y_hi_num

C ================================ Executable Code =============================

C Initialize some variables

!	write (6,*) x_fp, y_fp, dx_fp, dy_fp, cerflag
	hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick

C Initialize ok_hut to zero

	ok_hut = .false.
c
	resmult= 1.0

! TH - DC resolution based on pionCT analysis for SOS (average value ~1.8)
! this takes into account the effect of multiple tracks at high rates.
!	resmult= 1.8

C Initialize the xdc and ydc arrays to zero

	do i=1,12
	  xdc(i) = 0.
	  ydc(i) = 0.
	enddo

C------------------------------------------------------------------------------C
C                           Top of loop through hut                            C
C------------------------------------------------------------------------------C

C Go to the Cerenkov. (Drift back from focal plane). Cherenkov has 0.5 atm 
C of Freon, and is located 4.0 meters behind the focal point.
c The Cherenkov is 3m long, and has an Aluminum entrance window.
c only have cerenkov when p > 7500 MeV
c
c	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
c set decay flag false to avoid double counting

        if (cer_flag) then
	drift = hcer_1_zentrance
	call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	   radw = hfoil_exit_thick/hfoil_exit_radlen
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	 radw = hcer_entr_thick/hcer_entr_radlen
	 if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	 drift = hcer_1_zmirror - hcer_1_zentrance-hcer_mirglass_thick/2
	 radw = drift/hcer_1_radlen
	call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	if(ms_flag ) call musc_ext(m2,p,radw,drift,
     > dydzs,dxdzs,ys,xs)


	radw = hcer_mirglass_thick/hcer_mirglass_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
c
	drift = hcer_1_zexit - hcer_1_zmirror -hcer_mirglass_thick/2
	radw = drift/hcer_1_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,
     > dydzs,dxdzs,ys,xs)
	radw = hcer_exit_thick/hcer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        
       else
c  if vaccum pipe drift to pipe exit which is at same zpos as the cerenkov exit window. 
          if (vac_flag) then
   	   drift =  hcer_1_zexit 
 	   call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	   radw = hfoil_exit_thick/hfoil_exit_radlen
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          else ! helium bag
   	   drift = hcer_1_zentrance 
 	   call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	   radw = hfoil_exit_thick/hfoil_exit_radlen
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   radw = helbag_al_thick/helbag_al_radlen ! assume no distance to Helium bag Al mylar
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   radw = helbag_mylar_thick/helbag_mylar_radlen ! assume no distance to Helium bag Al mylar
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   drift = hcer_1_zexit - hcer_1_zentrance
	   radw = drift/helbag_hel_radlen
           call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	   if(ms_flag) call musc_ext(m2,p,radw,drift,
     > dydzs,dxdzs,ys,xs)
	   radw = helbag_al_thick/helbag_al_radlen ! assume no distance to Helium bag Al mylar
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	   radw = helbag_mylar_thick/helbag_mylar_radlen ! assume no distance to Helium bag Al mylar
	   if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          endif
        endif


C Go to first drift chamber set
C For simplicity, perform air MS (probably negligeable) at before drift
C instead of 1/2 way through.

*	write (6,*) 'made it to the first drift chamber'
	drift = (hdc_1_zpos - 0.5*hdc_nr_plan*hdc_del_plane) 
     > - hcer_1_zexit
	radw = drift/hair_radlen
	call project(xs,ys,drift,.false.,dflag,m2,p,pathlen)
	if (ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	jchamber = 1
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane) = xs + hdc_sigma(npl_off+iplane)
     > *tmpran1*resmult
	  ydc(npl_off+iplane) = ys + hdc_sigma(npl_off+iplane)
     > *tmpran2*resmult
	  if (iplane.eq.2 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
c	     write(*,*) ' iplane = ',iplane
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif
*	  write (6,*) 'i am still in the first drift chamber'
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	if (xs.gt.(hdc_1_bot-hdc_1x_offset) .or.
     >      xs.lt.(hdc_1_top-hdc_1x_offset) .or.
     >      ys.gt.(hdc_1_left-hdc_1y_offset) .or.
     >      ys.lt.(hdc_1_right-hdc_1y_offset) ) then
	  shmsSTOP_dc1 = shmsSTOP_dc1 + 1
c	  write(6,*) 'Lost in DC! delta,xp,yp',dpps,dxdzs,dydzs
	  if (use_det_cut) then
	     shmsSTOP_id = 34
	     goto 500
	  endif
	endif
        spec(39)=xs
        spec(40)=ys
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of first drift chamber set, drift to next

*	write (6,*) 'made it to the second drift chamber in the hut'
	drift = hdc_2_zpos - hdc_1_zpos - hdc_nr_plan*hdc_del_plane
C       Break this into 2 parts to properly account for decay.
C       We've already done decay up to the half-way point between the chambers.
	call project(xs,ys,drift/2.0,.false.,dflag,m2,p,pathlen)
	call project(xs,ys,drift/2.0,decay_flag,dflag,m2,p,pathlen)
	radw = drift/hair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	jchamber = 2
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane) = xs + hdc_sigma(npl_off+iplane)
     > *tmpran1*resmult
	  ydc(npl_off+iplane) = ys + hdc_sigma(npl_off+iplane)
     > *tmpran2*resmult
	  if (iplane.eq.2 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	if (xs.gt.(hdc_2_bot-hdc_2x_offset) .or.
     >      xs.lt.(hdc_2_top-hdc_2x_offset) .or.
     >      ys.gt.(hdc_2_left-hdc_2y_offset) .or.
     >      ys.lt.(hdc_2_right-hdc_2y_offset) ) then
	  shmsSTOP_dc2 = shmsSTOP_dc2 + 1
	  if (use_det_cut) then
	     shmsSTOP_id = 35
	     goto 500
	  endif
	endif
        spec(41)=xs
        spec(42)=ys
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of second drift chamber set, drift to the 1st hodoscope

	drift = hscin_1x_zpos - hdc_2_zpos - 0.5*hdc_nr_plan
     > *hdc_del_plane
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	xs_pad1x_x=xs !store the x position at pad1x                                                                                                                                                                                                                           
        xs_pad1x_y=ys
        xs_lo =  pscin_1x_center+pscin_1x_spacing*(pad_1x_lo_num-1)
     >      -pscin_1x_size/2.
        xs_hi =  pscin_1x_center+pscin_1x_spacing*(pad_1x_hi_num-1)
     >      +pscin_1x_size/2.
        hit_s1x = 0

	  if ( ys.le.(hscin_1x_left+hscin_1y_offset) .and.
     >      ys.ge.(hscin_1x_right+hscin_1y_offset)  .and.
     >      xs.le. xs_hi .and. xs.ge. xs_lo) then
            hit_s1x = 1
                 endif
 
c	if (ys.gt.(hscin_1x_left+hscin_1y_offset) .or.
c     >      ys.lt.(hscin_1x_right+hscin_1y_offset)) then
c	  shmsSTOP_s1 = shmsSTOP_s1 + 1
c	  if (use_det_cut) then
c	     shmsSTOP_id=36
c	     goto 500
c	  endif
c	endif

        spec(44)=ys
	radw = hscin_1x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	drift = hscin_1y_zpos - hscin_1x_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

 	pscin_1y_spacing = 7.5
        pscin_1y_center= +45.0 !-45.0                                                                                                                                                                                                                                          
        pscin_1y_size = 8.0 !8.0                                                                                                                                                                                                                                               
        ys_pad1y_y=ys !store y pos at 1y                                                                                                                                                                                                                                       
        ys_pad1y_x=xs
        ys_hi =  pscin_1y_center-pscin_1y_spacing*(pad_1y_lo_num-1)                                                                                                                                               
     >      +pscin_1y_size/2.
        ys_lo =  pscin_1y_center-pscin_1y_spacing*(pad_1y_hi_num-1)                                                                                                                                                
     >      -pscin_1y_size/2.
        hit_s1y = 0

        if (xs.le.(hscin_1y_bot+hscin_1x_offset) .and. ! was .or.!was xs.le                                                                                                                                                                                                    
     >      xs.ge.(hscin_1y_top+hscin_1x_offset) .and.
 	    ys .le. (ys_hi+0) .and. ys .ge. (ys_lo-0)) then
           	hit_s1y=1
		endif

     
c if (xs.gt.(hscin_1y_bot+hscin_1x_offset) .or.
c     >      xs.lt.(hscin_1y_top+hscin_1x_offset)) then
c	  shmsSTOP_s1 = shmsSTOP_s1 + 1
c	  if (use_det_cut) then
c	     shmsSTOP_id=37
c	     goto 500
c	  endif
c	endif
        spec(43)=xs
 	radw = hscin_1y_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to the second cherenkov

*	write (6,*) 'made it to the second cherenkov in the hut'
	drift = hcer_2_zentrance - hscin_1y_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = hcer_2_entr_thick/hcer_2_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_2_zmirror - hcer_2_zentrance
	radw = drift/hcer_2_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = hcer_mir_thick/hcer_mir_radlen

	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_2_zexit - hcer_2_zmirror
	radw = drift/hcer_2_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = hcer_2_exit_thick/hcer_2_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to 2nd hodoscope

	drift = hscin_2x_zpos - hcer_2_zexit
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

  	xs_pad2x_x=xs
        xs_pad2x_y=ys
        pscin_2x_center= -61.75
        pscin_2x_size = 10.0
        pscin_2x_spacing = 9.5
        xs_lo =  pscin_2x_center+pscin_2x_spacing*(pad_2x_lo_num-1)
     >      -pscin_2x_size/2.
        xs_hi =  pscin_2x_center+pscin_2x_spacing*(pad_2x_hi_num-1)
     >      +pscin_2x_size/2.
        hit_s2x=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                                                                                                                                                       
        if (ys.le.(hscin_2x_left+hscin_2y_offset) .and.                                                                                                                                                                                                    
     >      ys.ge.(hscin_2x_right+hscin_2y_offset) .and.
     >      xs.le.xs_hi .and. xs.ge .xs_lo) then                                                                                                                                                                                                                    
           hit_s2x=1
           endif

c	if (ys.gt.(hscin_2x_left+hscin_2y_offset) .or.
c     >      ys.lt.(hscin_2x_right+hscin_2y_offset)) then
c	  shmsSTOP_s3 = shmsSTOP_s3 + 1
c	  if (use_det_cut) then
c	     shmsSTOP_id = 38
c	     goto 500
c	  endif
c	endif
        spec(46)=ys
	radw = hscin_2x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	drift = hscin_2y_zpos - hscin_2x_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	ys_pad2y_y=ys
        ys_pad2y_x=xs
        pscin_2y_spacing = 5.0 !PARAM files                                                                                                                                                                                                                                    
        pscin_2y_center= 50.0 !-50.0                                                                                                                                                                                                                                                                                                                                                                                
        ys_hi =  pscin_2y_center-pscin_2y_spacing*(pad_2y_lo_num-1)!ys_lo=pscin_2y_center+pscin_2y_spacing*(pad_2y_lo_num-1)                                                                                                                                                   
     >      +pscin_2y_size/2.
        ys_lo =  pscin_2y_center-pscin_2y_spacing*(pad_2y_hi_num-1)!ys_hi=pscin_2y_center+pscin_2y_spacing*(pad_2y_hi_num-1)                                                                                                                                                   
     >      -pscin_2y_size_hi/2.

        hit_s2y=0
        if (xs.le.(hscin_2y_bot+hscin_2x_offset) .and. ! was .or.                                                                                                                                                                                                              
     >      xs.ge.(hscin_2y_top+hscin_2x_offset) .and. ! was then                                                                                                                                                                                                              
     >      ys .le. (ys_hi) .and. ys .ge. (ys_lo)) then    ! added Mike 5/10                                                                                                                                                                                                   
               hit_s2y=1
               endif

         if (use_det_cut
     >   .and.((hit_s1x+hit_s1y+hit_s2x+hit_s2y).lt.3)) then                                                                                                                                                                                                                
                shmsSTOP_id=39
                goto 500
             endif

 
c	if (xs.gt.(hscin_2y_bot+hscin_2x_offset) .or.
c     >      xs.lt.(hscin_2y_top+hscin_2x_offset)) then
c	  shmsSTOP_s2 = shmsSTOP_s2 + 1
c	  if (use_det_cut) then
c	     shmsSTOP_id=39
c	     goto 500
c	  endif
c	endif
        spec(45)=xs
	radw = hscin_2y_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Don't need to drift to calorimeter unless it's required in your trigger.
C Note that even with the standard PID trigger, the calorimeter is NOT
C required, since the trigger needs either the cerenkov OR the calorimeter.
C There is a seperate fiducial cut needed if you require the calorimeter
C in you analysis.  That cut is applied AFTER fitting the track (see below).

	drift = hcal_4ta_zpos - hscin_2y_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	if (ys.gt.hcal_left .or. ys.lt.hcal_right .or.
     >	   xs.gt.hcal_bottom .or. xs.lt.hcal_top) then
	  shmsSTOP_cal = shmsSTOP_cal + 1
	  if (use_det_cut) then
	     shmsSTOP_id=40
	     goto 500
	  endif
	endif
        spec(47)=xs
        spec(48)=ys




C and fit track to give new focal plane values, use LFIT from GENLIB

	do jchamber=1,hdc_nr_cham
	  npl_off = (jchamber-1)*hdc_nr_plan
	  do iplane=1,hdc_nr_plan
	    if (jchamber.eq.1) zdc(npl_off+iplane) = hdc_1_zpos +
     >          (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	    if (jchamber.eq.2) zdc(npl_off+iplane) = hdc_2_zpos +
     >          (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	  enddo
	enddo
	call lfit(zdc,xdc,12,0,dx,x0,badf)
	call lfit(zdc,ydc,12,0,dy,y0,badf)


	x_fp = x0
	y_fp = y0
	dx_fp = dx
	dy_fp = dy

C If you use a calorimeter cut in your analysis, the engine applied a
C a fiducial cut at the calorimeter.  This is based on the position of the
C TRACK at the calorimeter, not the real position of the event.  Go to the
C back of the calorimeter since engine uses a FID cut at the back.
C The standard fiducial cut is 5 cm from the edges of the block.

*	write (6,*) 'at the shower counter in the hut'
	xcal = x_fp + dx_fp * hcal_4ta_zpos
	ycal = y_fp + dy_fp * hcal_4ta_zpos
	if (ycal.gt.(hcal_left-5.0) .or. ycal.lt.(hcal_right+5.0) .or.
     >	   xcal.gt.(hcal_bottom-5.0) .or. xcal.lt.(hcal_top+5.0)) then
	  shmsSTOP_cal = shmsSTOP_cal + 1
	  if (use_det_cut) then
	     shmsSTOP_id=41
	     goto 500
	  endif
	endif
        spec(49)=xs
        spec(50)=ys

	ok_hut = .true.

C We are done with this event, whether GOOD or BAD.

500	continue

C ALL done!
	close (unit=chan)
	return
	end
