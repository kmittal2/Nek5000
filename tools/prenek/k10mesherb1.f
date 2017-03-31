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
      ITEM(4)='STASH CHUNK'
      NCHOIC = 4
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
      ELSEIF(CHOICE.EQ.'QUAD REGION')THEN
        call quadmesher
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

      call prs('SEE GUIDE - ENTER POINTS CCW$')
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
      subroutine makemesh(xmf,ymf,bcs,e1,e2)
      include 'basics.inc'
      integer e1,e2
      real xmf(e2,e1),ymf(e2,e1)
      integer i,e,j,f,ifld,dum
      character*3 bcs(4)

      nelold = nel
      write(6,*) nelold,'old element count'
         
      do j=1,e1-1
      do i=1,e2-1
        nel = nel+1
        x(1,nel) = xmf(i,j) 
        x(2,nel) = xmf(i+1,j) 
        x(3,nel) = xmf(i+1,j+1)  
        x(4,nel) = xmf(i,j+1) 
        y(1,nel) = ymf(i,j) 
        y(2,nel) = ymf(i+1,j) 
        y(3,nel) = ymf(i+1,j+1)  
        y(4,nel) = ymf(i,j+1) 
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
       write(6,*) f,e,bcs(4),'k10fe'
       cbc(f,e,ifld) = bcs(4)
      enddo
      enddo

      f=4
      do ifld=1,maxfld
      do e=nelold+1,nel,e2
       write(6,*) f,e,bcs(1),'k10fe'
       cbc(f,e,ifld) = bcs(1)
      enddo
      enddo

      f=3
      do ifld=1,maxfld
      do e=nelold+e1*(e2-1)-1,nel
       write(6,*) f,e,bcs(2),'k10fe'
       cbc(f,e,ifld) = bcs(2)
      enddo
      enddo

      f=2
      do ifld=1,maxfld
      do e=nelold+e2,nel,e2
       write(6,*) f,e,bcs(3),'k10fe'
       cbc(f,e,ifld) = bcs(3)
      enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine geomspace(a,b,n,r,xd)
      real a,b,r
      integer n,i,j
      real xd(100),dx,sc
 
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
