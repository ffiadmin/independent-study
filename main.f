c     Memory allocation parameters
      integer maxbas, maxnuc
      parameter (maxstp=500)
      parameter (maxbas=100)
      parameter (maxbs3=300)
      parameter (maxnuc=20)
      
c     Assign the "nbasis" variable to "static" memory
      integer nbasis
      common /basis_static/ nbasis
      
c     Assign the "nbasis" coordinate variables to "static" memory
      integer l(maxbas), m(maxbas), n(maxbas)
      double precision alpha(maxbas), x(maxbas), y(maxbas), z(maxbas)
      common /basis_coord_static/ alpha, l, m, n, x, y, z
      
c     Assign the "nnuc" variable to "static" memory
      integer nnuc
      common /nucleus_static/ nnuc
      
c     Assign the "nnuc" coordinate variable to "static" memory
      double precision nchg(maxnuc), nx(maxnuc), ny(maxnuc), nz(maxnuc)
      common /nucleus_coord_static/ nchg, nx, ny, nz
      
c     Assign additional program parameters to "static" memory
      integer nalp, nthet
      double precision alpstp, alpstrt, max, min, thstart, thstep
      common /param_static/ alpstp, alpstrt, max, min, nalp, nthet,
     +thstart, thstep
     
cc
c     Initialization
c     ---------------------------------
c
     
      call setNBasis(5)

c     In real life, call this function in a loop
c     do 10 i = 1, nbasis
c          call setBasisCoords(i, 1, 2, 3, 4.0D0, 5.0D0, 6.0D0, 7.0D0)
c10   continue

      call setBasisCoords(1, 0, 0, 0, 0.0D0, 0.0D0, 0.0D0,
     +0.190236205409D2)
      call setBasisCoords(2, 0, 0, 0, 0.0D0, 0.0D0, 0.0D0,
     +0.100124318636D2)
      call setBasisCoords(3, 0, 0, 0, 0.0D0, 0.0D0, 0.0D0,
     +0.526970098087D1)
      call setBasisCoords(4, 0, 0, 0, 0.0D0, 0.0D0, 0.0D0,
     +0.277352683204D1)
      call setBasisCoords(5, 0, 0, 0, 0.0D0, 0.0D0, 0.0D0,
     +0.145975096423D1)
      
      call setNNuc(3)

c     In real life, call this function in a loop
c     do 20 i = 1,nnuc
c          call setNucleusCoords(i, 1.0D0, 2.0D0, 3.0D0, 4.0D0)
c20   continue

      call setNucleusCoords(1, 0.0D0, 0.0D0, 0.0D0, 2.0D0)
      call setNucleusCoords(2, 0.0D0, 0.0D0, 0.4120D0, -1.0D0)
      call setNucleusCoords(3, 0.0D0, 0.0D0, -0.4120D0, -1.0D0)
      
      call setFinalParams(1, 1, 1.000D0, 0.00D0, 0.000D0, 0.00D0, 
     + -100, 1000000)
     
      call go()
      stop
      end
      
cc
c     Function library, section by Oliver Spryn
c     ------------------------------------------------------------------
c
      
cc
c     Set the value of the "nbasis" variable. The "nbasis"
c     variable is responsible for configuring the number of
c     basis functions.
c     
c     PRECONDITION:  C++ has opened an input text file.
c     POSTCONDITION: The "nbasis" variable will have been set.
c
c     REPLACES IN ORIGINAL PROGRAM:
c
c         read(5,*)nbasis
c         write(6,*) 'there are  ',nbasis,'  basis functions'
c
c     @access public
c     @param  integer value The value to assign to "nbasis"
c     @return void
c     @since  1.0.0
c

      subroutine setNBasis(value)
c     Import variables from "static" memory
          integer nbasis
          common /basis_static/ nbasis
          
c     Assign the subroutine parameter a type
          integer value
          
c     Assign "nbasis" its value
          nbasis = value
          
c     Alert the user of the program's progress
          print *, "There are", nbasis, " basis functions."
      end
      
cc
c     Set the coordinates for each of the nbasis values. The
c     "nbasis" value must have already been set. This is
c     designed to assign each of these values a coordinate.
c
c     The index (idx parameter) is NOT ZERO BASED!!! To access
c     the first element in a Fortran array, use array(1). Keep
c     that in mind when assigning this parameter a value. Thus,
c     the values passed into this parameter will range from:
c     1 ... nbasis, inclusive.
c
c     PRECONDITION:  setNBasis() has been called.
c     POSTCONDITION: Each nbasis value will have been assigned
c                    a set of coordinates
c
c     REPLACES IN ORIGINAL PROGRAM:
c
c         read(5,*)l(i),m(i),n(i),x(i),y(i),z(i),alpha(i)
c         x(i) = x(i)/0.529177249
c         y(i) = y(i)/0.529177249
c         z(i) = z(i)/0.529177249
c         write(6,*)l(i),m(i),n(i),x(i),y(i),z(i),alpha(i)
c     
c     @access public
c     @param  integer          idx      The index of the nbasis to assign these coordinates
c     @param  integer          lVal
c     @param  integer          mVal
c     @param  integer          nVal
c     @param  double precision xVal
c     @param  double precision yVal
c     @param  double precision zVal
c     @param  double precision alphaVal
c     @return void
c     @since  1.0.0
c

      subroutine setBasisCoords(idx, lVal, mVal, nVal, xVal,
     +yVal, zVal, alphaVal)
c     Import variables and arrays from "static" memory
          integer maxbas
          parameter (maxbas=100)
          
          integer l(maxbas), m(maxbas), n(maxbas)
          double precision alpha(maxbas), x(maxbas), y(maxbas),
     +z(maxbas)
          common /basis_coord_static/ alpha, l, m, n, x, y, z
          
c     Assign the subroutine parameters a type
          integer idx, lVal, mVal, nVal
          double precision alphaVal, xVal, yVal, zVal
      
c     Assign the nbasis index its respective values
          alpha(idx) = alphaVal
          l(idx) = lVal
          m(idx) = mVal
          n(idx) = nVal
          x(idx) = xVal
          y(idx) = yVal
          z(idx) = zVal
          
c     Convert the coordinates from angstroms to au
          x(idx) = x(idx)/0.529177249
          y(idx) = y(idx)/0.529177249
          z(idx) = z(idx)/0.529177249
          
c     Alert the user of the program's progress
          print *, "Assigned the coordinates for nbasis value:", idx
      end
      
cc
c     Set the value of the "nnuc" variable. The "nnuc"
c     variable is responsible for configuring the number of
c     nuclei.
c
c     PRECONDITION:  setBasisCoords() has been called.
c     POSTCONDITION: The "nnuc" variable will have been set.
c
c     REPLACES IN ORIGINAL PROGRAM:
c
c         read(5,*)nnuc
c         write(6,*)'there are  ',nnuc,'  nuclei.'
c     
c     @access public
c     @param  integer value The value to assign to "nnuc"
c     @return void
c     @since  1.0.0
c

      subroutine setNNuc(value)
c     Import variables from "static" memory
          integer nnuc
          common /nucleus_static/ nnuc
          
c     Assign the subroutine parameter a type
          integer value
          
c     Assign "nnuc" its value
          nnuc = value
          
c     Alert the user of the program's progress
          print *, "There are", nnuc, " nuclei."
      end

cc
c     Set the coordinates for each of the nnuc values. The
c     "nnuc" value must have already been set. This is
c     designed to assign each of these values a coordinate.
c
c     The index (idx parameter) is NOT ZERO BASED!!! To access
c     the first element in a Fortran array, use array(1). Keep
c     that in mind when assigning this parameter a value. Thus,
c     the values passed into this parameter will range from:
c     1 ... nnuc, inclusive.
c
c     PRECONDITION:  setNNuc() has been called.
c     POSTCONDITION: Each nnuc value will have been assigned
c                    a set of coordinates
c
c     REPLACES IN ORIGINAL PROGRAM:
c
c         read(5,*)nx(i),ny(i),nz(i),nchg(i)
c         nx(i) = nx(i)/0.529177249
c         ny(i) = ny(i)/0.529177249
c         nz(i) = nz(i)/0.529177249
c         write(6,*)nx(i),ny(i),nz(i),nchg(i)
c     
c     @access public
c     @param  integer          idx      The index of the nnuc to assign these coordinates
c     @param  double precision nChgVal
c     @param  double precision nxVal
c     @param  double precision nyVal
c     @param  double precision nzVal
c     @return void
c     @since  1.0.0
c

      subroutine setNucleusCoords(idx, nxVal, nyVal, nzVal, nChgVal)
c     Import variables and arrays from "static" memory
          integer maxnuc
          parameter (maxnuc=20)
          
          double precision nchg(maxnuc), nx(maxnuc), ny(maxnuc),
     +nz(maxnuc)
          common /nucleus_coord_static/ nchg, nx, ny, nz
          
c     Assign the subroutine parameters a type
          integer idx
          double precision nChgVal, nxVal, nyVal, nzVal
      
c     Assign the nnuc index its respective values
          nchg(idx) = nChgVal
          nx(idx) = nxVal
          ny(idx) = nyVal
          nz(idx) = nzVal
          
c     Convert the coordinates from angstroms to au
          nx(idx) = nx(idx)/0.529177249
          ny(idx) = ny(idx)/0.529177249
          nz(idx) = nz(idx)/0.529177249
          
c     Alert the user of the program's progress
          print *, "Assigned the coordinates for nnuc value:", idx
      end
      
cc
c     Set the value of the several final parameters to 
c     configure this program to begin processing.
c
c     PRECONDITION:  setNucleusCoords() has been called.
c     POSTCONDITION: A final set of configuration parameters
c                    will have been assigned a value.
c
c     REPLACES IN ORIGINAL PROGRAM:
c     
c         read(5,*)nthet,nalp,alpstrt,alpstp,thstart,thstep
c         ... several write statements ...
c         read(5,*)min,max
c         ... several write statements ...
c     
c     @access public
c     @param  integer          nThetVal
c     @param  integer          nAlpVal
c     @param  double precision alpStrtVal
c     @param  double precision alpStpVal
c     @param  double precision thStartVal
c     @param  double precision thStepVal
c     @param  integer          minVal
c     @param  integer          maxVal
c     @return void
c     @since  1.0.0
c

      subroutine setFinalParams(nThetVal, nAlpVal, alpStrtVal, 
     +alpStpVal, thStartVal, thStepVal, minVal, maxVal)
c     Import variables from "static" memory
          integer nalp, nthet
          double precision alpstp, alpstrt, max, min, thstart, thstep
          common /param_static/ alpstp, alpstrt, max, min, nalp, nthet,
     +thstart, thstep
          
c     Assign the subroutine parameters a type
          integer nAlpVal, nThetVal
          double precision alpStpVal, alpStrtVal, maxVal, minVal,
     +thStartVal, thStepVal
          
c     Assign the program parameters their respective values
          alpstp = alpStpVal
          alpstrt = alpStrtVal
          max = maxVal
          min = minVal
          nalp = nAlpVal
          nthet = nThetVal
          thstart = thStartVal
          thstep = thStepVal
          
c     Alert the user of the program's progress
          print *, "The initial values for r and theta are:", alpstrt,
     +"and " , thstart
          print *, "The program will take", nthet, " steps."
          print *, "Each step in theta will be", thstep, " radians."
          print *, "The program will take", nalp, " steps."
          print *, "Each step in r will be", alpstp, " au."
          print *, "The range for the real part of the energy is from",
     +min, "to ", max, "."
      end
      
cc
c     Written by: Dr. Michael Falcetta
c
c     Begin processing the data...
c
c     PRECONDITION:  setFinalParams() has been called.
c     POSTCONDITION: The program will begin its calculations.
c     
c     @access public
c     @return void
c     @since  1.0.0
c
      
      subroutine go()
c     Import variables and arrays from "static" memory
      integer maxbas, maxnuc
      parameter (maxstp=500)
      parameter (maxbas=100)
      parameter (maxbs3=300)
      parameter (maxnuc=20)
      
      integer nbasis
      common /basis_static/ nbasis
      
      integer l(maxbas), m(maxbas), n(maxbas)
      double precision alpha(maxbas), x(maxbas), y(maxbas), z(maxbas)
      common /basis_coord_static/ alpha, l, m, n, x, y, z
      
      integer nnuc
      common /nucleus_static/ nnuc
      
      double precision nchg(maxnuc), nx(maxnuc), ny(maxnuc), nz(maxnuc)
      common /nucleus_coord_static/ nchg, nx, ny, nz
      
      integer max, min, nalp, nthet
      double precision alpstp, alpstrt, thstart, thstep
      common /param_static/ alpstp, alpstrt, max, min, nalp, nthet,
     +thstart, thstep
      
c     Continue with the calculations...
      double precision PI,t1,t2,t3,t4,norm(maxbas),fact2
      double precision ovrlp,r1
      double precision kenergy

      double precision smat(maxbas,maxbas),
     + rr1,rr2,alpreal,theta,svec(maxbas),
     + aux(maxbs3),rwork(maxbs3),
     + dermag,dermin,ang(maxstp),angmin

      double complex scale,vint(maxbas,maxbas),
     + tint(maxbas,maxbas),xmatx(maxbas,maxbas),
     + xmat(maxbas,maxbas),core(maxbas,maxbas),
     + xtmp(maxbas,maxbas),XX1,XX2,eigval(maxbas),
     + VL(maxbas,maxbas),VR(maxbas,maxbas),
     + work(maxbs3),venergy,c6,eta(maxstp),
     + val(maxstp),deriv(maxstp),etamin,valmin

      integer naux,info,kk

      PI=dacos(-1.0d+00)

      scale = 1.0d+00
      XX1=(1.0d+00,0.0d+00)
      XX2=(0.0d+00,0.0d+00)
      
c     Calculate normalization constants

      do 30 i = 1,nbasis
           t1 = 2**(l(i)+m(i)+n(i))
           t2 = alpha(i)**(2*l(i)+2*m(i)+2*n(i)+3)
           t2 = t2**0.25
           t3 = (2.0d+00/PI)**0.75
           t4 = fact2(2*l(i)-1)*fact2(2*m(i)-1)*fact2(2*n(i)-1)
           t4= dsqrt(t4)
           norm(i) = t1*t2*t3/t4
30    continue

      do 11 ir=1,nbasis
           do 12 jc=1,nbasis
                call overlap(alpha(ir),l(ir)
     +          ,m(ir),n(ir),x(ir),y(ir),z(ir),
     +          alpha(jc),l(jc),m(jc),n(jc),
     +          x(jc),y(jc),z(jc),ovrlp)

                r1 = norm(ir)*norm(jc)*ovrlp
                smat(ir,jc)=r1                   


                call kinetic
     +          (alpha(ir),l(ir),m(ir),n(ir),x(ir),y(ir),z(ir),
     +          alpha(jc),l(jc),m(jc),n(jc),x(jc),y(jc),z(jc),
     +           kenergy)

                r1 =norm(ir)*norm(jc)*kenergy
                tint(ir,jc)=dcmplx(r1,0.0d+00)
12         continue
11     continue

      info=0
      naux=maxbs3
      call dsyev('V','U',nbasis,smat,maxbas,svec,aux,naux,info)

      do 98 iii = 1,nbasis
c     write(6,*)svec(iii)
98    continue

      alpreal=alpstrt

      do 888 istep=1,nalp 

      theta = thstart

      do 99 nscan=1,nthet
           rr1=alpreal*cos(theta)
           rr2=-alpreal*sin(theta)

           scale = dcmplx(rr1,rr2)
c          write(6,*)'scale ',scale


      do 40 ir=1,nbasis
           do 50 jc=1,nbasis


                vint(ir,jc)=(0.0d+00,0.0d+00)

              do 60 kk=1,nnuc     
                    call attrct(
     +  x(ir),y(ir),z(ir),norm(ir),l(ir),m(ir),n(ir),alpha(ir),
     +  x(jc),y(jc),z(jc),norm(jc),l(jc),m(jc),n(jc),alpha(jc),
     +  nx(kk),ny(kk),nz(kk),venergy,scale)

                    c6 = dcmplx(nchg(kk),0.0d+00)
                    vint(ir,jc)=vint(ir,jc)+c6*venergy
c       write(6,*) 'row, col, V ',ir,jc,vint(ir,jc),venergy,nchg(kk)

60              continue



                core(ir,jc)=scale*scale*tint(ir,jc)+vint(ir,jc)

c               write(6,*) 'row, col, S ',ir,jc,smat(ir,jc)
c               write(6,*) 'row, col, T ',ir,jc,tint(ir,jc)
c               write(6,*)'core  ',ir,jc,core(ir,jc)

50         continue
40    continue

C     Having built the ie matricies find the eigenvalues/vects of Smat


      do 70 ir = 1,nbasis
          do 80 jc = 1,nbasis
              r1 = smat(ir,jc)/dsqrt(svec(jc))
              xmat(ir,jc)=dcmplx(r1,0.0d+00)
              xmatx(jc,ir)=xmat(ir,jc)
c         write(6,*)xmat(ir,jc),xmatx(ir,jc)
80        continue
70    continue

      call zgemm('N','N',nbasis,nbasis,nbasis,XX1,core,maxbas,
     +   xmat,maxbas,XX2,xtmp,maxbas)

      call zgemm('N','N',nbasis,nbasis,nbasis,XX1,xmatx,maxbas,
     +   xtmp,maxbas,XX2,core,maxbas)

      call zgeev('V','V',nbasis,core,maxbas,eigval,VL,maxbas,
     +   VR,maxbas,work,maxbs3,rwork,info)

      do 86 iii = 1,nbasis
            eigval(iii)=(27.2114d+00,0.0d+00)*eigval(iii)
            ctest=dreal(eigval(iii))
            if ((ctest.ge.min).and.(ctest.le.max)) then
c           write(6,*)alpreal,theta,scale,eigval(iii)
            eta(nscan)=scale
            val(nscan)=eigval(iii)
            ang(nscan)=theta
            endif
86    continue

c     call zgemm('N','N',nbasis,nbasis,nbasis,XX1,VR,maxbas,
c    +   xmatx,maxbas,XX2,xtmp,maxbas)

c     call zgemm('N','N',nbasis,nbasis,nbasis,XX1,xmat,maxbas,
c    +   xtmp,maxbas,XX2,VR,maxbas)

c     do 88 jc =1,nbasis
c        do 89 ir = 1,nbasis
c           write(6,*)ir,jc,VR(ir,jc)
c89       continue
c88    continue
      
      
      theta=theta+thstep

99    continue

      dermin=1.0d+00

      do 999 i=1,nthet-1
          deriv(i)=(val(i+1)-val(i))/
     +    (eta(i+1)-eta(i))
          deriv(i)=deriv(i)
          dermag=dreal(deriv(i)*dconjg(deriv(i)))
          if (dermag.lt.dermin) then
              dermin=dermag
              etamin=eta(i)
              valmin=val(i)
              angmin=ang(i)
          endif
999   continue

      write(6,*)alpreal,angmin,etamin,valmin,dermin

      alpreal=alpreal+alpstp

888   continue          
      stop
      end

cc
c     Function library, original by Dr. Michael Falcetta
c     ------------------------------------------------------------------
c

      double precision FUNCTION fact2(n)
C     double factorial function
      INTEGER n


      fact2=1.0d+00

      if ((n.eq.0).or.(n.eq.1)) fact2=1.0d+00
      if (n.eq.3) fact2=3.0d+00
      if (n.eq.5) fact2=15.0d+00
      if (n.eq.7) fact2=15.0d+00*7.0d+00
      if (n.eq.9) fact2=15.0d+00**7.0d+00*9.0d+00
c     write(6,*)n,fact2

      END

      SUBROUTINE overlap(a1,l1,m1,n1,x1,y1,z1,
     +    a2,l2,m2,n2,x2,y2,z2,ovrlp)

      double precision a1,a2,ovrlp,x1,x2,y1,y2,z1,z2
      integer l1,l2,m1,m2,n1,n2

      double precision rab2,gamma,xp,yp,zp,pre,wx,wy,wz,
     +PI,OL1D
    

      PI=dacos(-1.0d+00)

      rab2 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
      gamma = a1+a2

      xp = (a1*x1+a2*x2)/(a1+a2)
      yp = (a1*y1+a2*y2)/(a1+a2)
      zp = (a1*z1+a2*z2)/(a1+a2)

      pre = (PI/gamma)**1.5d+00*dexp(-a1*a2*rab2/gamma)

c     write(6,*)'xp,yp,zp,pre'
c     write(6,*)xp,yp,zp,pre

      wx = OL1D(l1,l2,xp-x1,xp-x2,gamma)
      wy = OL1D(m1,m2,yp-y1,yp-y2,gamma)
      wz = OL1D(n1,n2,zp-z1,zp-z2,gamma)

c     write(6,*)'wx,wy,wz'     
c     write(6,*)l1,l2,wx,wy,wz     

      ovrlp = pre*wx*wy*wz

      return
      end

      double precision FUNCTION OL1D(l1,l2,PAx,PAb,gamma)

      integer l1,l2,max
      double precision PAx, PAb, gamma, sum,fact2,
     +PBINM

c     write(6,*)'l1,l2,PAx,PAb,gamma'
c     write(6,*)l1,l2,PAx,PAb,gamma
      max = int(0.5*(l1+l2))
c     write(6,*)'max  ',max, l1,l2
      sum=0.0d+00

      do 10, i = 0,max

           sum = sum + PBINM(2*i,l1,l2,PAx,PAb)*
     +     fact2(2*i-1)/(2.0d+00*gamma)**i

10    continue
c          write(6,*)'sum ' ,sum

      OL1D = sum
      end

      double precision FUNCTION PBINM(s,ia,ib,xpa,xpb)    


      integer s,ia,ib
      double precision xpa,xpb,sum,BIN

c     write(6,*)'xpa,xpb',xpa,xpb
      sum = 0.0d+00

      do 10 i =0,s
            if(((s-ia).le.i).and.(i.le.ib)) then
                sum = sum + BIN(ia,s-i)*BIN(ib,i)*
     +          xpa**(ia-s+i)*xpb**(ib-i)
            endif

10    continue

      PBINM = sum
c     write(6,*)'PBINM',PBINM

      end

      double precision FUNCTION  BIN(a,b)

      integer a,b,c,t1,t2,t3,fct

      t1 = fct(a)
      t2=fct(b)
      t3=fct(a-b)

c     write(6,*)'t1,t2,t3,',t1,t2,t3

      c = fct(a)/(fct(b)*fct(a-b))

c     write(6,*)'a,b,c ',a,b,c
      BIN=dfloat(c)
c     write(6,*)'BIN**',a,b,BIN

      end

      integer FUNCTION fct(a)

      integer a,sum

      sum =1

      do 10  i=0,a
          if (i.ge.1)sum = sum*i
10     continue

c     write(6,*)'a,sum',a,sum
      fct = sum

      end


      SUBROUTINE kinetic
     + (a1,l1,m1,n1,x1,y1,z1,a2,l2,m2,n2,x2,y2,z2,kenergy)

      double precision a1,x1,y1,z1,a2,x2,y2,z2,kenergy,
     +                 OL1,OL2,OL3,OL4

      double precision term0,term1,term2

      integer l1,m1,n1,l2,m2,n2

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2,x2,y2,z2,OL1)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2+2,m2,n2,x2,y2,z2,OL2)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2+2,n2,x2,y2,z2,OL3)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2+2,x2,y2,z2,OL4)

      term0 = a2*OL1*dfloat(2*l2+2*m2+2*n2+3)

      term1= -2.0d+00*a2*a2*(OL2+OL3+OL4) 

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2-2,m2,n2,x2,y2,z2,OL2)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2-2,n2,x2,y2,z2,OL3)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2-2,x2,y2,z2,OL4)

      term2= -0.5*(dfloat(l2*(l2-1))*OL2 +
     +             dfloat(m2*(m2-1))*OL3 +
     +             dfloat(n2*(n2-1))*OL4 )

      kenergy=term0+term1+term2

      return
      end


      SUBROUTINE attrct(
     +           x1,y1,z1,norm1,l1,m1,n1,a1,                             
     +           x2,y2,z2,norm2,l2,m2,n2,a2,                             
     +           nx,ny,nz,venergy,scale)

      double precision x1,y1,z1,norm1,a1
      double precision x2,y2,z2,norm2,a2
      double precision nx,ny,nz
      double precision gamma,xp,yp,zp,rab2
      double precision Ax(20),Ay(20),Az(20) 
      double precision r1,r2,r3,PI,r7

      integer l1,m1,n1,l2,m2,n2

      double complex venergy,sum,rcp2,c7

      double complex c1,c2,c3,gtest,c4,scale,
     + cx,cy,cz

      PI=dacos(-1.0d+00)

      gamma=a1+a2

      xp = (a1*x1+a2*x2)/(a1+a2)
      yp = (a1*y1+a2*y2)/(a1+a2)
      zp = (a1*z1+a2*z2)/(a1+a2)

      cx = dcmplx(nx,0.0d+00)
      cy = dcmplx(ny,0.0d+00)
      cz = dcmplx(nz,0.0d+00)

      cx =dcmplx(xp,0.0d+00)-scale*nx 
      cy =dcmplx(yp,0.0d+00)-scale*ny 
      cz =dcmplx(zp,0.0d+00)-scale*nz 

      rab2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
c     rcp2 = cx*dconjg(cx) + cy*dconjg(cy) + cz*dconjg(cz)
      rcp2 = cx*cx + cy*cy + cz*cz

      call Aarray(l1,l2,xp-x1,xp-x2,xp-nx,gamma,Ax)
      call Aarray(m1,m2,yp-y1,yp-y2,yp-ny,gamma,Ay)
      call Aarray(n1,n2,zp-z1,zp-z2,zp-nz,gamma,Az)

c      do 555 i = 1,20
c     write(6,*)'Ax,Ay,Az',Ax(i),Ay(i),Az(i)
c555   continue

      sum=(0.0d+00,0.0d+00)
      c4=(0.0d+00,0.0d+00)

      do 10 i = 0,l1+l2
           do 20 j=0,m1+m2
                do 30 k = 0,n1+n2
                     r1 = Ax(i+1)*Ay(j+1)*Az(k+1)
                     r2=dfloat(i+j+k)+0.5d+00
                     c1=dcmplx(r2,0.0d+00)
                     c2 = rcp2*dcmplx(gamma,0.0d+00)
                     r3 = dsqrt(dreal(c2 * dconjg(c2)))
                     if(dabs(r3).lt.(1.0d-17))c3=(1.0d-17,0.0)
                     c3=gtest(r2,c2)        
                     sum =sum+dcmplx(r1,0.0d+00)*c3 
c     write(6,*)Ax(i+1),Ay(j+1),Az(k+1)
c     write(6,*)'i,j,k,sum ',i,j,k,sum,c3,r4

30              continue
20         continue
10    continue


      r7 = -norm1*norm2*2.0d+00*PI/gamma
      r7 = dexp(-a1*a2*rab2/gamma)*r7
      c7 = dcmplx(r7,0.0d+00)
c     write(6,*)'venergy ',venergy
      venergy = scale*c7*sum                           

c     write(6,*)'venergy ',venergy

      return
      end

      SUBROUTINE Aarray(l1,l2,PA,PB,CP,g,A)

      double precision PA,PB,CP,g,A(20),t1
      integer l1,l2,Imax,lim1,lim2,II

c     write(6,*)'entering Aarray, l1,l2,PA,PB,CP,g'
c     write(6,*)l1,l2,PA,PB,CP,g

      Imax=l1+l2

      Do 10, i=1,20
          A(i)=0.0d+00
10    continue

      do 20 i=0,Imax
          lim1 = int(dfloat(i/2))

          do 30 j = 0,lim1
              lim2=int(dfloat(i-2*j)/2.0d+00)

              do 40 k = 0,lim2
                  II = i-2*j-k+1
                  call Aterm(i,j,k,l1,l2,PA,PB,CP,g,t1)
c                 write(6,*)'return Aterm i,r,u,l1,l2,t1,II'
c                 write(6,*)i,j,k,l1,l2,t1,II
                  A(II)=A(II)+t1
c                 write(6,*)'A(I),I  ',A(II),II
40            continue
30        continue
20    continue
c     write(6,*)'leaving  Aarray, l1,l2,PA,PB,CP,g,A'
c     write(6,*)l1,l2,PA,PB,CP,g
c     do 133 kk=1,5
c         write(6,*)A(kk)
c133   continue

      return
      end

      SUBROUTINE Aterm(i,j,k,l1,l2,PA,PB,CP,g,t1)

      integer i,j,k,l1,l2,fct
      double precision PA,PB,CP,g,t1,PBINM

      t1 = (-1)**i
      t1 = t1*PBINM(i,l1,l2,PA,PB)
      t1 = t1*(-1)**k
      t1 = t1 * fct(i)
      t1 = t1 * CP**(i-2*j-2*k)
      t1 = t1 * (0.25/g)**(j+k)
      t1 = t1 / (fct(j)*fct(k)*fct(i-2*j-2*k))

c     write(6,*)'i,r,u,t1',i,j,k,t1

      return
      end

      double Complex FUNCTION gtest(a,x)
      double complex x, sum,c1,c2,c3,sum2,xx2,cdig,b
      double precision a,r1,r2,sumnrm
      integer n,q,fct

      sum = (0.0d+00,0.0d+00)
      xx2 = (0.0d+00,0.0d+00)
      sumnrm = dsqrt(dreal(x * dconjg(x)))
c     write(6,*)'sumnrm',sumnrm 
      if(sumnrm.lt.1.0d+00)then
          do 10 n=0,20
              r1=a+dfloat(n)
              c1 = dcmplx(r1,0.0d+00)
              q = fct(n)
              r2=dfloat(q)
              c2=dcmplx(r2,0.0d+00)
              c3=-x
              sum = sum + (c3)**n/(c1*c2)
10        continue
          sum=  (0.5d+00,0.0d+00)*sum
c         write(6,*)'int',sum   
          gtest=sum
          return
      else
c         write(6,*)a,x
          b=dcmplx(a,0.0d+00)
          sum = cdig(b,x)
c         write(6,*)'sum',sum
          sum2 = cdig(b,xx2)
c         write(6,*)'sum2',sum2
          sum = (0.5d+00,0.0d+00)*x**(-a)*(sum2-sum)
c         write(6,*)'cdig',sum   
          gtest = sum
          return
      endif
      return
      end

      double complex function cdig(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex alpha,x,cdh
      double complex re,one,p,q
      double precision xlim,zero,dnrm
      data re,one/0.36787944117144232,1./
      data xlim,zero/1.,0./ibuf/34/
c --- If x is near the negative real axis, then shift to x=1.
      if(dnrm(x).lt.xlim.or.dreal(x).lt.zero.and.
     + dabs(dimag(x)).lt.xlim)then
      cdig=re/cdh(alpha,one)
      ilim=dreal(x/re)
      do 1 i=0,ibuf-ilim
      call term(alpha,x,i,p,q)
      cdig=cdig+p*q
1     continue
      else

      cdig=cdexp(-x+alpha*cdlog(x))/cdh(alpha,x)
      endif
      return
      end
c
      subroutine term(alpha,x,i,p,q)
c --- Calculate p*q = -1**i(1-x**(alpha+i))/(alpha+i)i! carefully.
      double complex alpha,x,p,q,ci,alphai
      double complex zero,one,two,cdlx
      double precision tol,xlim,dnrm
      data zero,one,two/0.,1.,2./tol/3.d-7/xlim/39./
      if(i.eq.0)q=one
      ci=i
      alphai=alpha+ci
      if(x.eq.zero)then
      p=one/alphai
      if(i.ne.0)q=-q/ci
      return
      endif
      cdlx=cdlog(x)
c --- If (1-x**alphai)=-x**alphai on the computer,
c --- then change the inductive scheme to avoid overflow.
      if(dreal(alphai*cdlx).gt.xlim.and.i.ne.0)then
      p=p*(alphai-one)/alphai
      q=-q*x/ci
      return
      endif
      if(dnrm(alphai).gt.tol)then
      p=(one-x**alphai)/alphai
      else
      p=-cdlx*(one+cdlx*alphai/two)
      endif
      if(i.ne.0)q=-q/ci
      return
      end
c
      double complex function cdh(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex alpha,x,cdhs
      double complex one,term,sum,cn,alpha1
      double precision buf
      data one/1./buf/0./
c --- If Re(alpha-x) is too big, shift alpha.
      n=dreal(alpha-x)-buf
      if(n.gt.0)then
      cn=n+1
      alpha1=alpha-cn
      term=one/x
      sum=term
      do 1 i=1,n
      cn=n-i+1
      term=term*(alpha1+cn)/x
      sum=term+sum

1     continue
      sum=sum+term*alpha1/cdhs(alpha1,x)
      cdh=one/sum
      else
      cdh=cdhs(alpha,x)
      endif
      return
      end
c
      double complex function cdhs(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex zero,half,one,alpha,x
      double complex p0,q0,p1,q1,r0,r1,ci,factor
      double precision tol1,tol2,error,dnrm
      data zero,half,one/0.,0.5,1./
      data tol1,tol2,error/1.d10,1.d-10,5.d-18/ilim/100000/
      q0=one
      q1=one
      p0=x
      p1=x+one-alpha
      do 1 i=1,ilim
      ci=i
      if(p0.ne.zero.and.q0.ne.zero.and.q1.ne.zero)then
      r0=p0/q0
      r1=p1/q1
      if(dnrm(r0-r1).le.dnrm(r1)*error)then
      cdhs=r1
      return
      endif
c --------- Occasionally renormalize the sequences to avoid over(under)flow.
      if(dnrm(p0).gt.tol1.or.dnrm(p0).lt.tol2.or.
     * dnrm(q0).gt.tol1.or.dnrm(q0).lt.tol2)then
      factor=p0*q0
      p0=p0/factor
      q0=q0/factor
      p1=p1/factor
      q1=q1/factor
      endif
      endif
      p0=x*p1+ci*p0
      q0=x*q1+ci*q0
      p1=p0+(ci+one-alpha)*p1
      q1=q0+(ci+one-alpha)*q1
1     continue
c --- If the peripheral routines are written correctly,
c --- the following four statements should never be executed.
      cdhs=half*(r0+r1)
      return
      end
c
      double precision function dnrm(z)

      double complex z
      dnrm=dabs(dreal(z))+dabs(dimag(z))
      return
      end
