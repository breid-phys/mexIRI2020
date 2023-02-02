#include "fintrf.h"

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
c	  MATLAB mex wrapper for the IRI 2016 model
c	  
c	  Accepts as input a size [9,n] array containing inputs for
c	  the NeQuick model. Each column should have the following:
c	  
C     h:     height (km)
C     alat:  gg. latitude  (degrees N)
C     along: gg. longitude (degrees E)
C     mth:   month (1 .. 12)
C     flx:   10.7 cm solar radio flux (flux units)
C     UT:    Universal Time (hours)
c
c	  Must be run from the folder containing the ccirxx.asc
c	  and modip.asc files
c	  
C Output is size (n,1) electron density in units of m^-3
C   (electrons per cubic meter)
c	  
c	  To compile MEX file:
c	  mex -compatibleArrayDims -DO3 -DfPIC -Dstatic ...
c		-output mexNeQuick2 mex_nequick.F NeQuick_2.F
c
c	  NeQuick_2.F is a slightly modified version of NeQuick_2.for
c	  the FORTRAN STOP call causes MATLAB to crash, so it must be removed
c	  
c	  Ben Reid 2020
c


C     Declarations
C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      
C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN
      integer*4 mexPrintf

C     Local Variables      
      mwPointer in_ptr, outf_ptr, oar_ptr, uin_ptr      
      real*8 ijf(50), input(9), uparam(16)
      real*4 outf(20,1000),oar(100)
      LOGICAL*4 jf(50)
      integer nummax, numhei, tnmax
      real*4 heibeg, heiend, heistp
      logical UINPUT

C     Array information:
      mwPointer mrows, ncols, umrows, uncols
      mwSize asize, oarsize,outfsize
      COMMON /ERRUNIT/ERRSTOP
	  
      INTEGER*8 i,j,k

      if(nrhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:IRI2016:nInput',
     &                           'Two inputs required.')
      elseif(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:IRI2016:nOutput',
     &                           'Too many output arguments.')
      endif


c     get inputs
      in_ptr = mxGetPr(prhs(1))
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,ijf,asize)

c     copy over jf array
      DO 2737 I=1,50
        if(ijf(i) .le. 0.0D0) then
	  jf(i) = .false.
	else
	  jf(i) = .true.
	endif
2737  continue

c     writes off
      jf(38)=.true.
c     messages off
      jf(34)=.false.

c     user inputs
      UINPUT = .false.
      if( nrhs .gt. 2) then
	UINPUT = .true.
	uin_ptr = mxGetPr(prhs(3))
        umrows = mxGetM(prhs(3))
        uncols = mxGetN(prhs(3))
	if(umrows.ne.size(uparam)) then
         call mexErrMsgIdAndTxt ('MATLAB:IRI2020:numrows',
     &    'User parameter array must have size [16,n].')
	endif
      elseif(.not.all(
c		Do not run if requested user inputs are not provided
     &  JF((/8,9,10,13,14,15,16,17,25,27,32,43,44,45,46/)))) then
	  call mexErrMsgIdAndTxt ('MATLAB:IRI2020:uimiss',
     &   'User Inputs requested but not provided')
      endif

      nout = 0
      if(jf(1)) nout=nout+1
      if(jf(2)) nout=nout+3
      if(jf(3)) nout=nout+5
      if(jf(3).and..not.jf(6)) nout=nout+2
      if(.not.jf(24)) nout=nout+1

      in_ptr = mxGetPr(prhs(2))
      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      
      if(mrows.ne.size(input)) then
        call mexErrMsgIdAndTxt ('MATLAB:IRI2020:nmrows',
     &   'Date/Location array must have size [9,n].')
      endif

      if(UINPUT.and.ncols.ne.uncols) then
        call mexErrMsgIdAndTxt ('MATLAB:IRI2020:uimismatch',
     &  'Date/Location array and user parameter array must match')
      endif

      asize = 1*mrows
      tnmax = 1
      nummax=1000
      
      DO 1167 I=1,ncols
        call mxCopyPtrToReal8(in_ptr+(8*mrows*(i-1)),input,asize)

	heibeg = sngl(input(7))
	heiend = sngl(input(8))
	heistp = sngl(input(9))
		
        numhei=int(abs(heiend-heibeg)/abs(heistp))+1
        if(tnmax.lt.numhei) tnmax=numhei
1167  continue

	if(tnmax.gt.nummax) then
          call mexErrMsgIdAndTxt ('MATLAB:IRI2020:nummax',
     &    'Too many elevations requested.')
	endif
	nummax = tnmax

c      prepare output arrays
       outfsize = size(outf);
c      plhs(1) = mxCreateDoubleMatrix(size(outf,1)*ncols,nummax,0)
       plhs(1) = mxCreateDoubleMatrix(nummax,nout*ncols,0)
       outf_ptr = mxGetPr(plhs(1))

       if(nlhs .gt. 1) then
         plhs(2) = mxCreateDoubleMatrix(size(oar),1*ncols,0)
	 oar_ptr = mxGetPr(plhs(2))
         oarsize=size(oar)
       endif

c     required by iri_sub
      call read_ig_rz
      call readapf107
	
c     LOOP OVER ALL INPUTS	
      DO 1 I=1,ncols

      DO 6249 j=1,100
6249    oar(j)=-1.0      

      if(UINPUT) then
	call mxCopyPtrToReal8(uin_ptr+(8*umrows*(i-1)),
     &	 uparam,size(uparam))
     
C    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
        if(.not.jf(8)) then
          oar(1)=sngl(uparam(1))
        endif

C    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
        if(.not.jf(9)) then
          oar(2)=sngl(uparam(2))
        endif

C    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
C       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not 
C       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
        if(.not.jf(10)) then
          oar(15)=sngl(uparam(3))
		  oar(16)=sngl(uparam(4))
        endif

C    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
        if(.not.jf(13)) then
          oar(3)=sngl(uparam(5))
        endif

C    jf(14) =.false.     OARR(4)=user input for hmF1/km
        if(.not.jf(14)) then
          oar(4)=sngl(uparam(6))
        endif

C    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
        if(.not.jf(15)) then
          oar(5)=sngl(uparam(7))
        endif

C    jf(16) =.false.     OARR(6)=user input for hmE/km
        if(.not.jf(16)) then
          oar(6)=sngl(uparam(8))
        endif

C    jf(17) =.false.     OARR(33)=user input for Rz12
        if(.not.jf(17)) then
          oar(33)=sngl(uparam(9))
        endif

C    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
        if(.not.jf(25)) then
          oar(41)=sngl(uparam(10))
        endif

C    jf(27) =.false.     OARR(39)=user input for IG12
        if(.not.jf(27)) then
          oar(39)=sngl(uparam(11))
        endif

C    jf(32) =.false.     OARR(46)=user input for 81-day avg F10.7
        if(.not.jf(32)) then
          oar(46)=sngl(uparam(12))
        endif

C    jf(43) =.false.     OARR(10)=user input for B0
        if(.not.jf(43)) then
          oar(10)=sngl(uparam(13))
        endif

C    jf(44) =.false.     OARR(35)=user input for B1
        if(.not.jf(44)) then
          oar(35)=sngl(uparam(14))
        endif

C    jf(45) =.false.     OARR(89)=user input for HNEA (Ne lower boundary)
        if(.not.jf(45)) then
          oar(89)=sngl(uparam(15))
        endif

C    jf(46) =.false.     OARR(90)=user input for HNEE (Ne upper boundary)
        if(.not.jf(46)) then
          oar(90)=sngl(uparam(16))
        endif

      endif

      call mxCopyPtrToReal8(in_ptr+(8*mrows*(i-1)),input,asize)

C     Clear error flag
      ERRSTOP = -1.0

      call IRI_SUB(jf,int(input(1)),sngl(input(2)),sngl(input(3)),
     &	   	int(input(4)),int(input(5)),sngl(input(6)),
     &		sngl(input(7)),
     &	   	sngl(input(8)),sngl(input(9)),outf,OAR)

      if(ERRSTOP.gt.0) then
	call mexErrMsgIdAndTxt ('MATLAB:IRI2020:StopCalled',
     &   'IRI_SUB called a stop command. Attempting to exit safely.')
      endif

      k=0
c	  output Ne if requested
C               OUTF(1,*)  ELECTRON DENSITY/M-3
      if(jf(1)) then
	call mxCopyReal8ToPtr(dble(outf(1,1:nummax)),
     &   outf_ptr+8*nummax*((i-1)*nout+k),nummax)
	k=k+1
      endif

c	  output Ti/Te if requested
C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
C               OUTF(3,*)  ION TEMPERATURE/K
C               OUTF(4,*)  ELECTRON TEMPERATURE/K
      if(jf(2)) then
        do 43 j=2,4
       	  call mxCopyReal8ToPtr(dble(outf(j,1:nummax)),
     &     outf_ptr+8*nummax*((i-1)*nout+k),nummax)
	  k=k+1
43	continue
      endif

c	  output Ni if requested
C               OUTF(5,*)  O+ ION DENSITY/% or /M-3 if jf(22)=f 
C               OUTF(6,*)  H+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(7,*)  HE+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(8,*)  O2+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(9,*)  NO+ ION DENSITY/% or /M-3 if jf(22)=f
      if(jf(3)) then
        do 44 j=5,9
       	  call mxCopyReal8ToPtr(dble(outf(j,1:nummax)),
     &     outf_ptr+8*nummax*((i-1)*nout+k),nummax)
	  k=k+1
44	continue
      endif

c	  output Ni if requested
C                 AND, IF JF(6)=.FALSE.:
C               OUTF(10,*)  CLUSTER IONS DEN/% or /M-3 if jf(22)=f
C               OUTF(11,*)  N+ ION DENSITY/% or /M-3 if jf(22)=f
      if(jf(3).and..not.jf(6)) then
        do 45 j=10,11
          call mxCopyReal8ToPtr(dble(outf(j,1:nummax)),
     &      outf_ptr+8*nummax*((i-1)*nout+k),nummax)
	  k=k+1
45	continue
      endif

c	  output special D-Region values
C  if(jf(24)    OUTF(14,1:11) standard IRI-Ne for 60,65,..,110km 
C     =.false.)        12:22) Friedrich (FIRI) model at these heights 
C                      23:33) standard Danilov (SW=0, WA=0) 
C                      34:44) for minor Stratospheric Warming (SW=0.5) 
C                      45:55) for major Stratospheric Warming (SW=1) 
C                      56:66) weak Winter Anomaly (WA=0.5) conditions
C                      67:77) strong Winter Anomaly (WA=1) conditions
      if(.not.jf(24)) then
	call mxCopyReal8ToPtr(dble(outf(14,1:nummax)),
     &  outf_ptr+8*nummax*((i-1)*nout+k),nummax)
	k=k+1
      endif


C     OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS  
      if(nlhs .gt. 1) then
	call mxCopyReal8ToPtr(dble(oar),
     &  oar_ptr+(8*(i-1)*oarsize),oarsize)
      endif

1     CONTINUE

      return
      end

