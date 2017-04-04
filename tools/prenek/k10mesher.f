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
c UPDATE 4 - April 2 - 6:06 pm - boundary mesher is now complete
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
      ITEM(5)='CIRCLE IN RECTANGLE'
      ITEM(6)='SMOOTH'
      ITEM(7)='STASH CHUNK'
      NCHOIC = 7
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
      ELSEIF(CHOICE.EQ.'CIRCLE IN RECTANGLE')THEN
        call circleinrectangle
      ENDIF

      GO TO 1

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine circleinrectangle
      include 'basics.inc'
c 
c This creates the mesh outside a circle in a rectangle
c INPUTS - center of circle, Radius of circle 
c 4 corners of Rectangle CW along with Elements and spacing 
c Circle surface is assumed to be a WALL
      real xycir(2),xyrec(4,2),r12(2),radcir
      integer i,j,k,l,e12(2)
      character*3 bcs(4)

      call prs('SEE GUIDE - CIRCLE IN A BOUNDARY$')
      CALL prs('Please enter the center of the circle:$')
      CALL rerr(xycir(1),xycir(2))
      CALL prs('Please enter the radius of the circle:$')
      CALL rer(radcir)
      CALL prs('Enter the corners of rectangle CW:$')
      CALL prs('Enter the (X,Y) FOR 1st point:$')
      CALL rerr(xyrec(1,1),xyrec(1,2))
      CALL prs('Enter the (X,Y) FOR 2st point:$')
      CALL rerr(xyrec(2,1),xyrec(2,2))
      CALL prs('Enter the (X,Y) FOR 3st point:$')
      CALL rerr(xyrec(3,1),xyrec(3,2))
      CALL prs('Enter the (X,Y) FOR 4st point:$')
      CALL rerr(xyrec(4,1),xyrec(4,2))
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

c We have all the inputs now, we are ready to mesh this :)

      




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
      real r12(2)
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

      call makemeshfromcorners(xx,yy,e12(2)+1,e12(1)+1,r12,bcs)

      return
      end
c-----------------------------------------------------------------------
      subroutine makemeshfromcorners(xx,yy,nrr,ncc,r12,bcs)
      integer nrr,ncc,i,j,k,l
      real mdxd(nrr,ncc),mdyd(nrr,ncc)
      real xydum(4,2),xx(4),yy(4)
      real sv1(ncc),sv2(nrr),r12(2)
      real s1(ncc),s3(ncc),s2(nrr),s4(nrr)
      real dx,dy,tdx,tdy,f1,f2
      character*3 bcs(4)

      do i=1,4
       xydum(i,1) = xx(i)
       xydum(i,2) = yy(i)
      enddo
      call genquad(xydum,ncc,nrr,mdxd,mdyd,r12)

      call makemesh(mdxd,mdyd,bcs,ncc,nrr)
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

      call genbndrelmat(e12,r12,xyspl,nnpts,bcs) !xdum,ydum are in common block

c      call prexit
      goto 400
  300 CONTINUE
      write(6,*) 'There was error in reading the dat file'
      write(6,*) 'file name you said was ',fnamef
      call prexit
  400 CONTINUE
 
      return
      end
c-----------------------------------------------------------------------
      subroutine genbndrelmat(e12,r12,xyspl,nnpts,bcs)
c INPUT e12, r12, xyspl - indices differentiation between different edges
      PARAMETER (NMAX=10000)
      common /ctmpk10/ xdum(nmax),ydum(nmax),xn(nmax),yn(nmax),sc(nmax)
      integer e1,e2,e12(2),xyspl(4),nnpts,n1,n2,ind1,i,j
      real x1(nnpts),y1(nnpts)
      real x2(nnpts),y2(nnpts)
      real x3(nnpts),y3(nnpts)
      real x4(nnpts),y4(nnpts)
      real r12(2)
      real mfx(e12(2)+1,e12(1)+1),mdx(e12(2)+1,e12(1)+1)
      real mdx2(e12(2)+1,e12(1)+1)
      real mfy(e12(2)+1,e12(1)+1),mdy(e12(2)+1,e12(1)+1)
      real mdy2(e12(2)+1,e12(1)+1)
      character bcs(4)*3
c  mfx will be used first for x and then for y
c  mdx will be used for dummy matrix with zeros

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
      call rzero(mdx,e1*e2)
      call rzero(mdy,e1*e2)
c  LAST ROW - A-B
      j = e2+1
      do i=1,e1+1
       mdx(j,i) = x1(i)
       mdy(j,i) = y1(i)
      enddo
c  LAST COLUMN C - B down
      j=e1+1
      do i=1,e2+1
       mdx(i,j) = x2(e2+2-i) 
       mdy(i,j) = y2(e2+2-i) 
      enddo
c   TOP ROW D - C
      j = 1
      do i=1,e1+1
        mdx(j,i) = x3(e1+2-i)
        mdy(j,i) = y3(e1+2-i)
      enddo
c   FIRST COLUMN D - A DOWN
      j=1
      do i=1,e2+1
        mdx(i,j) = x4(i)
        mdy(i,j) = y4(i)
      enddo
c     MDX AND MDY ARE READY NOW 

c      do i=1,e2+1
c      do j=1,e1+1
c       write(6,*) i,j,mdx(i,j),mdy(i,j),'k10ij'
c      enddo
c      enddo

c     NOW NEED TO INTERPOLATE THIS TO INTERIOR
      call makeinteriormesh(mdx,mdy,e2+1,e1+1,mdx2,mdy2,r12)

      call changemeshcw(mdx2,mdy2,e2+1,e1+1)

      call makemesh(mdx2,mdy2,bcs,e1+1,e2+1)
      call redraw_mesh


      
      return
      end
c-----------------------------------------------------------------------
      subroutine changemeshcw(mdx2,mdy2,nrr,ncc)
      integer nrr,ncc,i,j,k,l
      real mdx2(nrr,ncc),mdy2(nrr,ncc),dumx(nrr,ncc)
c Takes a mesh which is input CCW and makes it CW
c ESSENTIALLY MAKES FIRST ROW -> LAST, 2nd row -> 2nd last etc...
  
      call copy(dumx,mdx2,nrr*ncc)
      do i=1,ncc
      do j=1,nrr
        mdx2(j,i) = dumx(nrr-j+1,i)
      enddo
      enddo

      call copy(dumx,mdy2,nrr*ncc)
      do i=1,ncc
      do j=1,nrr
        mdy2(j,i) = dumx(nrr-j+1,i)
      enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine makeinteriormesh(mdx,mdy,nrr,ncc,mdx2,mdy2,r12)
      integer nrr,ncc,i,j,k,l
      real mdx(nrr,ncc),mdy(nrr,ncc)
      real mdx2(nrr,ncc),mdy2(nrr,ncc)
      real mdxd(nrr,ncc),mdyd(nrr,ncc)
      real xydum(4,2),xydum1(ncc,2),xydum2(nrr,2)
      real sv1(ncc),sv2(nrr),r12(2)
      real s1(ncc),s3(ncc),s2(nrr),s4(nrr)
      real dx,dy,tdx,tdy,f1,f2
c this takes the boundary mesh in mdum and interpolates it to interior
c in md2.
      xydum(1,1) = mdx(1,1)
      xydum(2,1) = mdx(1,ncc)
      xydum(3,1) = mdx(nrr,ncc)
      xydum(4,1) = mdx(nrr,1)
      xydum(1,2) = mdy(1,1)
      xydum(2,2) = mdy(1,ncc)
      xydum(3,2) = mdy(nrr,ncc)
      xydum(4,2) = mdy(nrr,1)

      call genquad(xydum,ncc,nrr,mdxd,mdyd,r12)
c This genereates the mesh based on vertices

      do i=1,ncc
        xydum1(i,1) = mdxd(1,i)
        xydum1(i,2) = mdyd(1,i)
      enddo
      call getlenscale(xydum1,s1,0,ncc)

      do i=1,ncc
        xydum1(i,1) = mdxd(nrr,i)
        xydum1(i,2) = mdyd(nrr,i)
      enddo
      call getlenscale(xydum1,s3,0,ncc)

      do i=1,nrr
        xydum2(i,1) = mdxd(i,ncc)
        xydum2(i,2) = mdyd(i,ncc)
      enddo
      call getlenscale(xydum2,s2,0,nrr)

      do i=1,nrr
        xydum2(i,1) = mdxd(i,1)
        xydum2(i,2) = mdyd(i,1)
      enddo
      call getlenscale(xydum2,s4,0,nrr)

      call copy(mdx2,mdxd,nrr*ncc) 
      call copy(mdy2,mdyd,nrr*ncc) 

c MOVE EDGE 1
      do j=1,ncc
       dx = mdx(1,j) - mdxd(1,j);
       dy = mdy(1,j) - mdyd(1,j);
       do i=1,nrr
         f1 = (s1(ncc)-s1(j))/s1(ncc)
         f2 = (s2(nrr)-s2(i))/s2(nrr);
         tdx = dx*f2;
         tdy = dy*f2;
         mdx2(i,j) = mdx2(i,j)+tdx;
         mdy2(i,j) = mdy2(i,j)+tdy;
       enddo
      enddo


c MOVE EDGE 2
      do i=1,nrr
        dx = mdx(i,ncc) - mdxd(i,ncc);
        dy = mdy(i,ncc) - mdyd(i,ncc);
        do j=1,ncc
          f1 = (s1(j))/s1(ncc);
          f2 = (s2(nrr)-s2(i))/s2(nrr);
          tdx = dx*f1;
          tdy = dy*f1;
          mdx2(i,j) = mdx2(i,j)+tdx;
          mdy2(i,j) = mdy2(i,j)+tdy;
        enddo
       enddo

c MOVE EDGE 3
      do j=1,ncc
       dx = mdx(nrr,j) - mdxd(nrr,j);
       dy = mdy(nrr,j) - mdyd(nrr,j);
       do i=1,nrr
         f1 = (s1(j))/s1(ncc)
         f2 = (s2(i))/s2(nrr);
         tdx = dx*f2;
         tdy = dy*f2;
         mdx2(i,j) = mdx2(i,j)+tdx;
         mdy2(i,j) = mdy2(i,j)+tdy;
       enddo
      enddo


c MOVE EDGE 4
      do i=1,nrr
        dx = mdx(i,1) - mdxd(i,1);
        dy = mdy(i,1) - mdyd(i,1);
        do j=1,ncc
          f1 = (s1(ncc)-s1(j))/s1(ncc);
          f2 = (s2(i))/s2(nrr);
          tdx = dx*f1;
          tdy = dy*f1;
          mdx2(i,j) = mdx2(i,j)+tdx;
          mdy2(i,j) = mdy2(i,j)+tdy;
        enddo
       enddo

 
      return
      end
c-----------------------------------------------------------------------
      subroutine genquad(xydum,ncc,nrr,mdxd,mdyd,r12)
      integer nrr,ncc,i,j,k,l
      real xd1(ncc),xd2(nrr),r12(2)
      real xydum(4,2),mdxd(nrr,ncc),mdyd(nrr,ncc)
      real jx1(ncc,2),jx1t(2,ncc),dum1(ncc*nrr)
      real jx2(nrr,2),jx2t(2,nrr)
      real xyintd(4)
      integer ptr(4)


      call geomspace(0.,1.,ncc,r12(1),xd1)
      call geomspace(0.,1.,nrr,r12(2),xd2)

      xyintd(1) = 0.
      xyintd(2) = 1.
      call gen_int_gz(jx1,jx1t,xd1,ncc,xyintd,2)
      call gen_int_gz(jx2,jx2t,xd2,nrr,xyintd,2)

      ptr(1) = 1
      ptr(2) = 4
      ptr(3) = 2
      ptr(4) = 3

      call swapvecinds(xyintd,xydum(1,1),ptr,4)
      call mxm(jx2,nrr,xyintd,2,dum1,2)
      call mxm(dum1,nrr,jx1t,2,mdxd,ncc)

      call swapvecinds(xyintd,xydum(1,2),ptr,4)
      call mxm(jx2,nrr,xyintd,2,dum1,2)
      call mxm(dum1,nrr,jx1t,2,mdyd,ncc)
       
      
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
      subroutine getlenscale(xyd,sv2,flag,n)
      integer flag,i,n
      real xyd(n,2),sv2(n)

      sv2(1) = 0
      do i=1,n-1
       sv2(i+1) = sv2(i) +
     $            dist2d(xyd(i,1),xyd(i,2),xyd(i+1,1),xyd(i+1,2))
      enddo

      if (flag.eq.1) then
       call rescale(sv2,0.,1.,n)
      endif


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
