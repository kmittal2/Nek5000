C                       Copyright (C) 1990, by the
C
C               Massachusetts Institute of Technology  and Nektonics, Inc.
C
C All Rights Reserved
C
C This program is a licenced product of MIT and Nektonics, Inc.,  and it is
C not to be disclosed to others, copied, distributed, or displayed
C without prior authorization.
C
C------------------------------------------------------------------------------
C
C
C THE PURPOSE OF THIS FILE IS TO CREATE MESHES IN A MODULAR WAY
c UPDATE 1 - MARCH 30 - 2:46 AM - quad mesher works for 4 sided regions
c UPDATE 2 - MARCH 31 - 1:00 AM - boundary mesher is in progress. I have added
c most of the machinery. Need to mimic genquadmesh from matlab now
c a lot of other supporting functions have been added
c UPDATE 3 - March 31 - 4:51 PM - found bug in last night's code.. fixed that
c and started adding the boundary mesher part
c 
C
C
C
C
      subroutine k10meshmenu
      include 'basics.inc'
      integer icalld,e
      save    icalld
      data    icalld /0/
      logical iftmp
C
      iftmp=ifgrid
      ifgrid=.false.

1     continue
      ITEM(1)='BUILD MENU'
      ITEM(2)='DELETE'
      ITEM(3)='QUAD REGION'
      ITEM(4)='QUAD BOUNDARY'
      ITEM(5)='SMOOTH'
      ITEM(6)='STASH CHUNK'
      NCHOIC = 6
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'K10 MESHER')
3013  CONTINUE

      IF(CHOICE.EQ.'BUILD MENU')THEN
            NCURVE=0
            DO 114 IE=1,NEL
            DO 114 IEDGE=1,12
               IF (CCURVE(IEDGE,IE).NE.' ') THEN
                  NCURVE=NCURVE+1
               ENDIF
114         CONTINUE
            IFGRID=IFTMP
            RETURN
      ELSEIF(CHOICE.EQ.'DELETE') THEN
        call delete
      ELSE IF(CHOICE.EQ.'SMOOTH')THEN
         call smoother2(nel)
      ELSEIF(CHOICE.EQ.'QUAD REGION')THEN
        call quadmesher
      ELSEIF(CHOICE.EQ.'QUAD BOUNDARY')THEN
        call quadfromboundary
      ENDIF

      GO TO 1

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine quadmesher
      include 'basics.inc'
c multiple ways to make a mesh
c input x,y pairs ccw and specify nelx and nely
c also specify spacing - uniform or geometric 
c also specify boundary conditions 
c 
      real xx(4),yy(4)
      integer e12(2),i,j
      real r12(2),xd(100),xd2(100)
      real jx1(200),jx1t(200),dum1(200)
      real jx2(200),jx2t(200),dum2(200)
      real xmf(10000),ymf(10000)
      real xyintd(4) 
      integer ptr(4)
      character*3 bcs(4)

      call prs('SEE GUIDE - ENTER POINTS CW$')
      CALL prs('Enter the (X,Y) FOR 1st point:$')
      CALL rerr(xx(1),yy(1))
      CALL prs('Enter the (X,Y) FOR 2st point:$')
      CALL rerr(xx(2),yy(2))
      CALL prs('Enter the (X,Y) FOR 3st point:$')
      CALL rerr(xx(3),yy(3))
      CALL prs('Enter the (X,Y) FOR 4st point:$')
      CALL rerr(xx(4),yy(4))
      call prs('Enter number of elements on 1st and 2nd side:$')
      call reii(e12(1),e12(2))
      call prs('Enter GEOMETRIC spacing for both sides:$')
      call rerr(r12(1),r12(2))
      call prs('Enter the BCs for side 1:$')
      call res(bcs(1),3)
      call prs('Enter the BCs for side 2:$')
      call res(bcs(2),3)
      call prs('Enter the BCs for side 3:$')
      call res(bcs(3),3)
      call prs('Enter the BCs for side 4:$')
      call res(bcs(4),3)
      write(6,*) bcs
      write(6,*) xx,yy,e12,r12,'k10xyer'

      call geomspace(0.,1.,e12(1)+1,r12(1),xd)
      call geomspace(0.,1.,e12(2)+1,r12(2),xd2)

      xyintd(1) = 0.
      xyintd(2) = 1.
      call gen_int_gz(jx1,jx1t,xd,e12(1)+1,xyintd,2)
      call gen_int_gz(jx2,jx2t,xd2,e12(2)+1,xyintd,2)

      ptr(1) = 1
      ptr(2) = 4
      ptr(3) = 2
      ptr(4) = 3

      call swapvecinds(xyintd,xx,ptr,4)      
      call mxm(jx2,e12(2)+1,xyintd,2,dum1,2)
      call mxm(dum1,e12(2)+1,jx1t,2,xmf,e12(1)+1)

      call swapvecinds(xyintd,yy,ptr,4)      
      call mxm(jx2,e12(2)+1,xyintd,2,dum1,2)
      call mxm(dum1,e12(2)+1,jx1t,2,ymf,e12(1)+1)

      call makemesh(xmf,ymf,bcs,(e12(1)+1),(e12(2)+1))
      call redraw_mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine quadfromboundary
      include 'basics.inc'
      PARAMETER (NMAX=10000)
      common /ctmpk10/ xdum(nmax),ydum(nmax),xn(nmax),yn(nmax),sc(nmax)
      real xdum,ydum,xn,yn,sc
      integer nnpts,xyspl(4),e12(2)
      real r12(2)
      character fnamef*8, string2*72, fname*4, bcs(4)*3
      integer n,i,n1,n2
c    This takes as input 1 big file which has all xy coordinates
c    the user must also tell what are the ranges for each side
c    example xy has 100 points so if the user says sides are between
c    1,20,60,100 means S1 is from 1-20, S2 is 20-60 and so on
c    sides must be in clocwise order
c    the user must also tell geometric spacing factor and e1 and e2
c    and also bcs
      call prs('SEE GUIDE - You need to specify file name with xy$')
      CALL prs('Enter the name of the comma separated dat file:$')
      CALL prs('The name can be only 4 characters long:$')
      CALL res(fnamef,4)
      fname = '.dat'
      fnamef(5:8) = fname

      CALL prs('Please enter the points seperating the 4 lines:$')
      call reiiii(xyspl(1),xyspl(2),xyspl(3),xyspl(4))

      call prs('Enter number of elements on 1st and 2nd side:$')
      call reii(e12(1),e12(2))
      call prs('Enter GEOMETRIC spacing for both sides:$')
      call rerr(r12(1),r12(2))
      call prs('Enter the BCs for side 1:$')
      call res(bcs(1),3)
      call prs('Enter the BCs for side 2:$')
      call res(bcs(2),3)
      call prs('Enter the BCs for side 3:$')
      call res(bcs(3),3)
      call prs('Enter the BCs for side 4:$')
      call res(bcs(4),3)

  
      write(6,*) fnamef,xyspl,'file name and critical points'
      open (unit=9,file=fnamef,status='old', iostat=ierr)

      do i=1,NMAX
        READ(9,*,ERR=300,END=200) xdum(i),ydum(i)
      enddo
  200 CONTINUE

      nnpts = i-1
      write(6,*) nnpts,' points have been read from dat file'
c     so far we have the data read. now we have to construct the matrix
c     etc...

      call genbndrelmat(e12,r12,xyspl,nnpts) !xdum,ydum are in common block

      call prexit

      n1 = 10
      n2 = 6
      write(6,*) 'try 1'
      call geomspace(0.,1.,n2,1.1,sc)
      call geomspace(3.,5.,n1,1.2,xdum) 
      call geomspace(1.,2.,n1,1.1,ydum) 

      call changescale(xn,yn,xdum,ydum,sc,n2,n1)

      do i=1,6
       write(6,*) i,sc(i),xn(i),yn(i),'k10ints'
      enddo
      do i=1,10
       write(6,*) i,xdum(i),ydum(i),'k10origs'
      enddo

      call prexit
      goto 400
  300 CONTINUE
      write(6,*) 'There was error in reading the dat file'
      write(6,*) 'file name you said was ',fnamef
      call prexit
  400 CONTINUE
 
      return
      end
c-----------------------------------------------------------------------
      subroutine genbndrelmat(e12,r12,xyspl,nnpts)
c INPUT e12, r12, xyspl - indices differentiation between different edges
      PARAMETER (NMAX=10000)
      common /ctmpk10/ xdum(nmax),ydum(nmax),xn(nmax),yn(nmax),sc(nmax)
      integer e1,e2,e12(2),xyspl(4),nnpts,n1,n2,ind1,i,j
      real x1(nnpts),y1(nnpts)
      real x2(nnpts),y2(nnpts)
      real x3(nnpts),y3(nnpts)
      real x4(nnpts),y4(nnpts)
      real r12(2)
      real mf(e12(2)+1,e12(1)+1),mdum(e12(2)+1,e12(1)+1)
c  mf will be used first for x and then for y
c  mdum will be used for dummy matrix with zeros

      e1 = e12(1)
      e2 = e12(2)

c    1st edge
      n1 = xyspl(2)-xyspl(1)+1
      n2 = e12(1)+1
      ind1 = xyspl(1)
      call copy(xn,xdum(ind1),n1)
      call copy(yn,ydum(ind1),n1)
      call geomspace(0.,1.,n2,r12(1),sc)
      call changescale(x1,y1,xn,yn,sc,n2,n1)
c    2nd edge
      n1 = xyspl(3)-xyspl(2)+1
      n2 = e12(2)+1
      ind1 = xyspl(2)
      call copy(xn,xdum(ind1),n1)
      call copy(yn,ydum(ind1),n1)
      call geomspace(0.,1.,n2,r12(2),sc)
      call changescale(x2,y2,xn,yn,sc,n2,n1)
c    3rd edge
      n1 = xyspl(4)-xyspl(3)+1
      n2 = e12(1)+1
      ind1 = xyspl(3)
      call copy(xn,xdum(ind1),n1)
      call copy(yn,ydum(ind1),n1)
      call geomspace(0.,1.,n2,r12(1),sc)
      call changescale(x3,y3,xn,yn,sc,n2,n1)
c    4th edge
      n1 = nnpts-xyspl(4)+1
      n2 = e12(2)+1
      ind1 = xyspl(4)
      if (abs(xdum(nnpts)-xdum(1)).lt.1e-6.and.
     $     abs(ydum(nnpts)-ydum(1)).lt.1e-6) then
c      that means the wrap is closed
      
      else
c      no periodic wrap..make it periodic
       n1=n1+1
       xdum(nnpts+1)=xdum(1)
       ydum(nnpts+1)=ydum(1)
      endif 
      call copy(xn,xdum(ind1),n1)
      call copy(yn,ydum(ind1),n1)
      call geomspace(0.,1.,n2,r12(2),sc)
      call changescale(x4,y4,xn,yn,sc,n2,n1)
c ALL EDGES HAVE BEEN DONE AT THIS POINT

c    Now need to make the interpolation matrix with zeros in interior
      call rzero(mdum,e1*e2)
c  LAST ROW - A-B
      j = e2+1
      do i=1,e1+1
       mdum(j,i) = x1(i)
      enddo
c  LAST COLUMN C - B down
      j=e1+1
      do i=1,e2+1
       mdum(i,j) = x2(e2+2-i) 
      enddo
c   TOP ROW D - C
      j = 1
      do i=1,e1+1
        mdum(j,i) = x3(e1+2-i)
      enddo
c   FIRST COLUMN D - A DOWN
      j=1
      do i=1,e2+1
        mdum(i,j) = x4(i)
      enddo
c     MDUM is ready now
c     NOW NEED TO INTERPOLATE THIS TO INTERIOR

      call prexit
       
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine changescale(xn,yn,xx,yy,sc,n1,n2)
c  this takes xx and yy and interpoaltes using cubic spline to xn,yn
c  xx,yy has n2 points, xn,yn has n1 points
      integer n1,n2
      real xx(n2),yy(n2),xn(n1),yn(n1),s(n2),sc(n1)
      real scb(n1)

      call copy(scb,sc,n1)

      s(1) = 0.
      do i=1,n2-1
       s(i+1) = s(i)+dist2d(xx(i),yy(i),xx(i+1),yy(i+1))
      enddo

      call rescale(sc,0.,s(n2),n1)
      call spfit(s,xx,sc,xn,n1,n2)
      call spfit(s,yy,sc,yn,n1,n2)

      call copy(sc,scb,n1)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine spfit(xx,yy,xn,yn,n1,n2)
c   given values of xin and yin, what are values of yout based on xout
c   n2 are points in xx, n1 are points in xn
      integer n,i,j,n1,n2
      real xx(n2),xn(n1),yy(n2),yn(n1),yp0,yp2(n2)

      yp0 = 0.0
      call spline(xx,yy,n2,yp0,yp0,yp2)

      do i=1,n1
          call splint(xx,yy,yp2,n2,xn(i),yn(i))
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine rescale(xx,x1,x2,n)
c rescales x from x1 to x2 
      integer i,n
      real xx(n),x1,x2,xn(n)
      real minv,maxv

      minv = glmin(xx,n)
      maxv = glmax(xx,n)
      do i=1,n
       xn(i) = (xx(i)-minv)/(maxv-minv)
      enddo

      do i=1,n
       xx(i) = x1+xn(i)*(x2-x1)
        if (i.eq.1) xx(i) = x1
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine makemesh(xmf,ymf,bcs,e1,e2)
      include 'basics.inc'
      integer e1,e2
      real xmf(e2,e1),ymf(e2,e1)
      real xmft(e1,e2),ymft(e1,e2)
      integer i,e,j,f,ifld,dum
      character*3 bcs(4)

      call transpose_r(xmft,e1,xmf,e2)
      call transpose_r(ymft,e1,ymf,e2)

      nelold = nel
      write(6,*) nelold,'old element count'
         
      do j=1,e2-1
      do i=1,e1-1
        nel = nel+1
        x(1,nel) = xmft(i,j) 
        x(2,nel) = xmft(i+1,j) 
        x(3,nel) = xmft(i+1,j+1)  
        x(4,nel) = xmft(i,j+1) 
        y(1,nel) = ymft(i,j) 
        y(2,nel) = ymft(i+1,j) 
        y(3,nel) = ymft(i+1,j+1)  
        y(4,nel) = ymft(i,j+1) 
      enddo
      enddo

      do ifld=1,maxfld
      do e=nelold+1,nel
      do f=1,4
        cbc(f,e,ifld) = 'E  '    ! totally arbitrary default
      enddo
      enddo
      enddo

      e1=e1-1
      e2=e2-1
      write(6,*) e1,e2,nel,'k10e1e2nel'
      f=1
      do ifld=1,maxfld
      do e=nelold+1,nelold+e2
       cbc(f,e,ifld) = bcs(1)
      enddo
      enddo

      f=2
      do ifld=1,maxfld
      do e=nelold+e1,nel,e1
       cbc(f,e,ifld) = bcs(2)
      enddo
      enddo

      f=3
      do ifld=1,maxfld
      do e=nelold+e1*(e2-1)+1,nel
       cbc(f,e,ifld) = bcs(3)
      enddo
      enddo

      f=4
      do ifld=1,maxfld
      do e=nelold+1,nel,e1
       cbc(f,e,ifld) = bcs(4)
      enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine geomspace(a,b,n,r,xd)
      real a,b,r
      integer n,i,j,norig
      real xd(n),dx,sc
c genereates n points between a and b with ratio r and stores in xd
      norig = n
      if (n.le.1) write(6,*) 'Invalid number of elements entered'
      if (n.le.1) call prexit
      if (b.ne.a) then
       n = n-1
       dx = (b-a)/n
       xd(1) = 0.
       do i=2,n+1
         xd(i) = xd(i-1)+dx
         dx = dx*r
       enddo 
 
       sc = (b-a)/(xd(n+1)-xd(1))
 
       do i=1,n+1
         xd(i) = a+sc*xd(i)
       enddo
       else
          do i=1,n+1
           xd(i) = a
         enddo
       endif
      n = norig
      return
      end
c-----------------------------------------------------------------------
      subroutine swapvecinds(vn,vo,inds,n)
      integer n,i
      integer inds(n)
      real vn(n),vo(n)
                
      do i=1,n
       vn(i) = vo(inds(i))
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C
      PARAMETER (NMAX=4000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
      Y2(1)=0.0
      U(1) =0.0
C
      DO 10 I=2,N-1
         IR=I+1
         IL=I-1
         SIG=(X(I)-X(IL))/(X(IR)-X(IL))
         P=SIG*Y2(IL)+2.
         Y2(I)=(SIG-1.)/P
         U(I)= ( 6.*
     $     ( (Y(IR)-Y(I))/(X(IR)-X(I))-(Y(I)-Y(IL))/ (X(I)-X(IL) ) )
     $            / (X(IR)-X(IL))
     $    - SIG*U(IL) )/P
   10 CONTINUE
C
      QN=0.0
      UN=0.0
C
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 20 K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
   20 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C     p. 88-89, numerical recipes
C
      DIMENSION XA(N),YA(N),Y2A(N)
C
      KLO=1
      KHI=N
    1   IF ((KHI-KLO).GT.1) THEN
           K=(KHI+KLO)/2
           IF (XA(K).GT.X) THEN
              KHI=K
           ELSE
              KLO=K
           ENDIF
           GOTO 1
        ENDIF
C
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0) THEN
         WRITE(6,*) XA(KHI), 'Hey buddy - you blew it.'
         RETURN
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
c-----------------------------------------------------------------------
