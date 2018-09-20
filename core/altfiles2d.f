      subroutine doaltschwarz_dummy2d(ngeomp,igeom)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'MVGEOM'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOP'
      INCLUDE 'CTIMER'
      INCLUDE 'GLOBALCOM'
      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELV)
     $ ,             RES2  (LX1,LY1,LZ1,LELV)
     $ ,             RES3  (LX1,LY1,LZ1,LELV)
     $ ,             DV1   (LX1,LY1,LZ1,LELV)
     $ ,             DV2   (LX1,LY1,LZ1,LELV)
     $ ,             DV3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)

      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP
      REAL           prcp   (LX2,LY2,LZ2,LELV)
      REAL           dprc   (LX2,LY2,LZ2,LELV)

      integer ngeomp,ntot
      integer idx1,idx2
      real             valint(lx1,ly1,lz1,lelt,nfldmax_nn)
      common /valmask/ valint
      common/ exprx/ exactdprx,exactdpry,exactpr
      real exactpr(lx2,ly2,lz2,lelv)
      real exactdprx(lx2,ly2,lz2,lelv)
      real exactdpry(lx2,ly2,lz2,lelv)
      real errpr(lx2,ly2,lz2,lelv)
      parameter (nptfpt = 100)
      real xlist(nptfpt),ylist(nptfpt),errint(nptfpt)
      real rwk(nptfpt,ldim+1)
      integer iwk(nptfpt,3)
      real wgt(lx1*ly1*lz1*lelt)
      common /globwgt / wgt

   
      ntot1 = lx1*ly1*lz1*nelv
      do i=1,ntot1
       xv = xm1(i,1,1,1)
       exactpr(i,1,1,1) = sin(x)*cos(y)
       exactdprx(i,1,1,1) = cos(x)*cos(y)
       exactdpry(i,1,1,1) = -sin(x)*sin(y)
      enddo
      call rzero(pr,lx1*ly1*lz1*nelv)
      rmeanexp = glsc2_univ(exactpr,wgt,ntot1)
     $          /glsc2_univ(wgt,wgt,ntot1)

c     Initial error
      call sub3(errpr,pr,exactpr,ntot1)
      call outpost(pr,errpr,vz,exactpr,pr,'   ')

      idx1 = 0
      idx2 = 1

      call modpresint('v  ','o  ')

      igeomps = 1
      if (nid.eq.0) write(6,*) idsess,ngeomp,' alt schwarz'
     
      do igeomp=1,ngeomp
           call neknek_xfer_fld(pr,ldim+1)
           call neknek_bcopy(ldim+1)
           call copy(prcp,pr,lx1*ly1*lz1*nelv)

ccc      Solve for session 1
           if (idsess.eq.idx1) then
           call crespsp_dummy2d(respr)
           call invers2  (h1,vtrans,ntot1)
           call rzero    (h2,ntot1)
           call ctolspl  (tolspl,respr)
           napproxp(1) = laxtp
           call hsolve   ('PRES',dpr,respr,h1,h2
     $                        ,pmask,vmult
     $                        ,imesh,tolspl,nmxh,1
     $                        ,approxp,napproxp,binvm1)
           call add2    (pr,dpr,ntot1)
           endif
           call neknekgsync()
ccc      Exchange data
           call neknek_xfer_fld(pr,ldim+1)
           call neknek_bcopy(ldim+1)
ccc      Solve for session 2
           if (idsess.eq.idx2) then
           call crespsp_dummy2d(respr)
           call invers2  (h1,vtrans,ntot1)
           call rzero    (h2,ntot1)
           call ctolspl  (tolspl,respr)
           napproxp(1) = laxtp
           call hsolve   ('PRES',dpr,respr,h1,h2
     $                        ,pmask,vmult
     $                        ,imesh,tolspl,nmxh,1
     $                        ,approxp,napproxp,binvm1)
           call add2    (pr,dpr,ntot1)
           endif
c           call ortho_univ2_tar   (pr,rmeanexp)
           call sub3(dprc,prcp,pr,ntot1)
           dprmax = uglamax(dprc,ntot1)

          call sub3(errpr,pr,exactpr,ntot1)
          call outpost(pr,errpr,vz,exactpr,pr,'   ')
          err_max = uglamax(errpr,ntot1)
         if (nid_global.eq.0)
     $      write(6,'(i2,i8,2i4,1p3e13.4,a11)') idsess,istep,igeom,
     $      igeomp,time,
     $      dprmax,err_max,' max-dp-nn'

         enddo
         call modpresint('o  ','v  ')

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine crespsp_dummy2d (respr)

C     Compute startresidual/right-hand-side in the pressure

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'GLOBALCOM'

      REAL           RESPR (LX1*LY1*LZ1,LELV)
c
      COMMON /SCRNS/ TA1   (LX1*LY1*LZ1,LELV)
     $ ,             TA2   (LX1*LY1*LZ1,LELV)
     $ ,             TA3   (LX1*LY1*LZ1,LELV)
     $ ,             WA1   (LX1*LY1*LZ1*LELV)
     $ ,             WA2   (LX1*LY1*LZ1*LELV)
     $ ,             WA3   (LX1*LY1*LZ1*LELV)
      COMMON /SCRMG/ W1    (LX1*LY1*LZ1,LELV)
     $ ,             W2    (LX1*LY1*LZ1,LELV)
     $ ,             W3    (LX1*LY1*LZ1,LELV)

      common /scruz/         sij (lx1*ly1*lz1,6,lelv)
      parameter (lr=lx1*ly1*lz1)
      common /scrvz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)
      common/ exprx/ exactdprx,exactdpry,exactpr
      real exactpr(lx2,ly2,lz2,lelv)
      real exactdprx(lx2,ly2,lz2,lelv)
      real exactdpry(lx2,ly2,lz2,lelv)


      CHARACTER CB*3

      NXYZ1  = lx1*ly1*lz1
      NTOT1  = NXYZ1*NELV
      NFACES = 2*ldim

      call rzero(W1,ntot1)
      do i=1,ntot1
       xv = xm1(i,1,1,1)
       yv = ym1(i,1,1,1)
       w1(i,1) = bm1(i,1,1,1)*2.*sin(xv)*cos(yv)
      enddo
      call dssum(w1)

c     add old pressure term because we solve for delta p 
      call rone(ta1,ntot1)
      call rzero   (ta2,ntot1)

      call bcdirpc (pr)
      call axhelm  (respr,pr,ta1,ta2,imesh,1)
c      call dssum(respr)

c      call col2(respr,binvm1,ntot1)
      call col2(w1,binvm1,ntot1)
      call outpost(respr,w1,vz,pr,t,'rhs')
      call col2(w1,bm1,ntot1)
c      call col2(respr,bm1,ntot1)

      call chsign  (respr,ntot1)
      call add2(respr,w1,ntot1)

      DO IEL=1,NELV
         DO IFC=1,2*ldim
            CALL RZERO  (W1(1,IEL),NXYZ1)
            CALL RZERO  (W2(1,IEL),NXYZ1)
            CB = CBC(IFC,IEL,IFIELD)
            IF (CB(1:1).EQ.'W'.OR.CB(1:1).EQ.'v') then
               CALL FACCL3
     $         (W1(1,IEL),exactdprx(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3
     $         (W2(1,IEL),exactdpry(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)

               CALL ADD2   (W1(1,IEL),W2(1,IEL),NXYZ1)
              IF (ldim.EQ.3)
     $         CALL ADD2   (W1(1,IEL),W3(1,IEL),NXYZ1)

               CALL FACCL2 (W1(1,IEL),AREA(1,1,IFC,IEL),IFC)
            ENDIF
            CALL ADD2 (RESPR(1,IEL),W1(1,IEL),NXYZ1)
       enddo
       enddo 

 

      return
      END
c----------------------------------------------------------------------
