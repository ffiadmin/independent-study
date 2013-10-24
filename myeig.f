      subroutine myeig(nbasis,nb2,Rmat,Imat,scale)


      integer nbasis,info,naux
      double precision Rmat(nbasis,nbasis),Imat(nbasis,nbasis),scale
      complex*16 Cval(nbasis),Cmat(nbasis,nbasis),t1,t2
      complex*16 VL(nbasis,nbasis),VR(nbasis,nbasis),WORK(100)
      double precision RWORK(nb2),aux(77),OUT(nbasis),t6,t7

      naux=77  
      do 10 irow=1,nbasis

            do 20 icol=1,nbasis
               
              t1 = dcmplx(0.0d+00,Imat(irow,icol))
              t1=scale*t1 
              t2 = dcmplx(Rmat(irow,icol),0.0d+00)
              Cmat(irow,icol)=t2-t1
C             write(6,*)irow,icol,Cmat(irow,icol)

20          continue
10    continue




C     call dsyev("V", "U", nbasis, Rmat, nbasis, OUT , aux, naux, info)
      call zgeev('N','V',nbasis,Cmat,nbasis,
     +   Cval,VL,nbasis,VR,nbasis,WORK,100,RWORK,INFO)

      write(6,*)"info =  ", info
      do 30 k = 1,nbasis
          t6=dimag(Cval(k))
          t7=dreal(Cval(k))
          t6=27.2114d+00*t6
          t7=27.2114d+00*t7

          write(6,*)scale, " ++ ", t7,t6    
C         write(6,*)"in myeig  ",scale, t7, t6 
30    continue

      return
      end
