      subroutine outfld2d_k10(u,v,w,p,tt,nlxy,name,ifld)
      include 'SIZE'
      include 'TOTAL'
      real u (lx1,ly1,lz1,lelt)
      real v (lx1,ly1,lz1,lelt)
      real w (lx1,ly1,lz1,lelt)
      real p (lx2,ly2,lz2,lelt)
      real tt(lx1,ly1,lz1,lelt)
      real work1(lx1,ly1)
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
     $    ,pa (lx1,ly2,lz2)     ,pb (lx1,ly1,lz2)
      common /scrns/ w1(lx1,ly1),w2(lx1,ly1),w3(lx1,ly1),w4(lx1,ly1)
     $ ,w5(lx1,ly1),w6(lx1,ly1)
      character*3 name
      integer eg,e,i,j,nlxy
      character*2  excode(15)
      character*20 outfile

      call mappr(pm1,p,pa,pb)

      call blank(excode,30)
c
      if(ifxyo) then
          excode(1) = 'X '
          excode(2) = 'Y '
      endif
      
      excode(3) = 'U '
      excode(4) = 'P '
      if(ifto) excode(5) = 'T '
      
      ierr = 0
      if (nid.eq.0) then
         call blank(outfile,20)
         if (ifld.lt.100) then
            write(outfile,2) name,ifld
    2       format(a3,'2d.fld',i2.2)
         elseif (ifld.lt.1000) then
            write(outfile,3) name,ifld
    3       format(a3,'2d.fld',i3)
         elseif (ifld.lt.10000) then
            write(outfile,4) name,ifld
    4       format(a3,'2d.fld',i4)
         elseif (ifld.lt.100000) then
            write(outfile,5) name,ifld
    5       format(a3,'2d.fld',i5)
         elseif (ifld.lt.1000000) then
            write(outfile,6) name,ifld
    6       format(a3,'2d.fld',i6)
         endif

        open(unit=24,file=outfile,status='unknown')
        call dump_header2d(excode,lx1,ly1,nlxy,1,ierr)
      endif


         do eg=1,nlxy
            call nekgsync()
            mid = gllnid(eg)
            call rzero(w1,lx1*ly1)
            call rzero(w2,lx1*ly1)
            call rzero(w3,lx1*ly1)
            call rzero(w4,lx1*ly1)
            call rzero(w5,lx1*ly1)
            call rzero(w6,lx1*ly1)
          if (mid.eq.nid) then
            e = gllel(eg)
             do j=1,ny1
             do i=1,nx1
             w1(i,j) = xm1(i,j,1,e)
             w2(i,j) = ym1(i,j,1,e)
             w3(i,j) = u  (i,j,1,e)
             w4(i,j) = v  (i,j,1,e)
             w5(i,j) = pm1(i,j,1,e)
             w6(i,j) = tt (i,j,1,e)
            enddo
            enddo
          endif
            call nekgsync()
            call gop(w1,work1,'+  ',lx1*ly1)
            call gop(w2,work1,'+  ',lx1*ly1)
            call gop(w3,work1,'+  ',lx1*ly1)
            call gop(w4,work1,'+  ',lx1*ly1)
            call gop(w5,work1,'+  ',lx1*ly1)
            call gop(w6,work1,'+  ',lx1*ly1)
            call nekgsync()
          if (nid.eq.0) then
           do j=1,ny1
           do i=1,nx1
          if((ifxyo).and.(ifto)) then
              write(24,'(1p6e14.6)') w1(i,j),w2(i,j),w3(i,j)
     $                          ,w4(i,j),w5(i,j),w6(i,j)
          elseif(ifxyo) then
               write(24,'(1p5e14.6)') w1(i,j),w2(i,j),w3(i,j)
     $                          ,w4(i,j),w5(i,j)
          elseif(ifto) then
                write(24,'(1p4e14.6)')w3(i,j)
     $                          ,w4(i,j),w5(i,j),w6(i,j) 
          else
                write(24,'(1p3e14.6)')w3(i,j)
     $                          ,w4(i,j),w5(i,j)
          endif 

           enddo
           enddo
          endif
          call nekgsync()
         enddo


         close(24)
      call err_chk(ierr,'Error using byte_write,outfld2d. Abort.$')
      return
      end
c-----------------------------------------------------------------------
      subroutine get_yavg_k10(xval)
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
     $    ,pa (lx1,ly2,lz2)     ,pb (lx1,ly1,lz2)
      real xval,dy2,diff,tol
      integer e,eg,i,j,k,n
      real vxav,w1,vyav,w2,pmav,w3
      real work1
      vxav = 0.
      vyav = 0.
      pmav = 0.
      w1 = 0.
      w2 = 0.
      w3 = 0.

      call mappr(pm1,pr,pa,pb)
      tol = 0.00005

      n = 0
      do e=1,nelt
       eg = lglel(e)
       k=1
       do i=1,nx1
         diff = abs(xm1(i,1,1,e)-xval)
         if (diff.lt.tol) then
         dy2 = .5*( ym1(i,ny1,k,e) - ym1(i,1,k,e) )
         do j=1,ny1
               n = n+1
               vxav = vxav+dy2*wym1(j)*vx(i,j,k,e)
               w1 = w1+dy2*wym1(j) ! redundant but clear
               vyav = vyav+dy2*wym1(j)*vy(i,j,k,e)
               w2 = w2+dy2*wym1(j) ! redundant but clear
               pmav = pmav+dy2*wym1(j)*pm1(i,j,k,e)
               w3 = w3+dy2*wym1(j) ! redundant but clear
         enddo
         endif
        enddo
      enddo

      call gop(vxav,work1,'+  ',1)
      call gop(w1,work1,'+  ',1)
      call gop(vyav,work1,'+  ',1)
      call gop(w2,work1,'+  ',1)
      call gop(pmav,work1,'+  ',1)
      call gop(w3,work1,'+  ',1)
      call gop(n,work1,'+  ',1)
      if (n.ne.0) then
        vxav = vxav/w1
        vyav = vyav/w2
        pmav = pmav/w3
        vmag = (vxav**2+vyav**2)**0.5
        if (nid.eq.0) then
      write(6,*) xval,vxav,vyav,vmag,pmav,'k10_xm1_vx_vy_vmag_pm'
        endif
      else
        if (nid.eq.0) then
          write(6,*) 'x-coordinate didnt match',xval
        endif
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine getpres_k10(nelxy)  ! This routine to modify mesh coordinates
      include 'SIZE'
      include 'TOTAL'

      integer z,e,f,t,ix,iy,iz,ex,ey,ez,eg
      parameter (lt=lx1*ly1*lz1*lelv)
      common /scrcg/ pm1(lx1,ly1,lz1,lelv),w1(lt), w2(lt)  ! Mapped pressure
      real cm1,cm2,cm3
      real dumz
      integer nelxy

      nelz  = nelgt/nelxy
      ntest = nelz*nelxy
      if (ntest.ne.nelgt) call exitti('nelxy fail:',nelxy)

      call mappr(pm1,pr,w1,w2)

      z = 0
      do eg=1,nelgt
         call nekgsync ! Puts all processors in lock-step
         mid=gllnid(eg)
         call nekgsync ! Puts all processors in lock-step
         if (mid.eq.nid) then
           call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
           e=gllel(eg)
           do f=1,4
             if (ez.eq.1.and.cbc(f,e,1).eq.'W  ') then
                call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
                iz = 1 ! Take 2nd plane
                do iy=ky1,ky2
                do ix=kx1,kx2
                  z=z+1
                  cm1 = xm1(ix,iy,iz,e)
                  cm2 = ym1(ix,iy,iz,e)
                  cm3 = pm1(ix,iy,iz,e)
                  write(6,2) cm1,cm2,cm3
                enddo
                enddo
              endif
            enddo
         endif
      enddo
    2 format(1p3e12.4,' xyp')
      return
      end
c-----------------------------------------------------------------------
      subroutine gettau_k10(nelxy)
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelv)
      common /scrns/ vort(lt,3),w1(lt),w2(lt)
      real vort1(lx1,ly1,lz1,lelv,3)
      integer z,e,f,t,ix,iy,iz,eg,ex,ey,ez
      real cm1,cm2,cm3,cm4,cm5
      real dumz,vis
      integer nelxy

      nelz  = nelgt/nelxy
      ntest = nelz*nelxy
      if (ntest.ne.nelgt) call exitti('nelxy fail:',nelxy)

      n = nx1*ny1*nz1*nelv
      call comp_vort3(vort1, w1, w2, vx, vy, vz)
      ifto = .true.
      ifxyo = .true.

      vis = param(2)

      do eg=1,nelgt
         call nekgsync ! Puts all processors in lock-step
         mid=gllnid(eg)
         call nekgsync ! Puts all processors in lock-step
         call nekgsync ! Puts all processors in lock-step
         if (mid.eq.nid) then
           call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
           e = gllel(eg)
           do f=1,4
             if (ez.eq.1.and.cbc(f,e,1).eq.'W  ') then
             call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
             iz = 2
             do iy=ky1,ky2
             do ix=kx1,kx2
               cm1 = xm1(ix,iy,iz,e)
               cm2 = ym1(ix,iy,iz,e)
               cm3 = vort1(ix,iy,iz,e,1)
               cm4 = vort1(ix,iy,iz,e,2)
               cm5 = vort1(ix,iy,iz,e,3)
               write(6,4) cm1,cm2,cm3,cm4,cm5
             enddo
             enddo
         endif
        enddo
       endif
      enddo
    4 format(1p5e12.4,' xyt')

      return
      end
c-----------------------------------------------------------------------
      subroutine z_avg(ua,u,gs_avg_hndl,nelxy,ifld)
      include 'SIZE'
      include 'TOTAL'

c     Compute the z average of quantity u() - assumes global tens.prod.


c     common usage
c      integer gs_avg_hndl
c      save    gs_avg_hndl
c      data    gs_avg_hndl / 0 /
c      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
c     $    ,pa (lx1,ly2,lz2)     ,pb (lx1,ly1,lz2)
c      real vxa(lx1*ly1*lz1*lelt)
c      real vya(lx1*ly1*lz1*lelt)
c      real vza(lx1*ly1*lz1*lelt)
c      real pma(lx1*ly1*lz1*lelt)
c
c
c     nelxy = 160
c     ifld  = 1
c      call mappr(pm1,p,pa,pb)
c      call z_avg(vxa,vx,gs_avg_hndl,nelxy,ifld)
c      call z_avg(vya,vy,gs_avg_hndl,nelxy,ifld)
c      call z_avg(vza,vz,gs_avg_hndl,nelxy,ifld)
c      call z_avg(pma,pm1,gs_avg_hndl,nelxy,ifld)

c      call outpost(vxa,vya,vza,pma,t,'   ')

      real u (lx1,ly1,lz1,lelt)
      real ua(lx1,ly1,lz1,lelt)

      integer gs_avg_hndl,e,ex,ey,ez,eg

      if (gs_avg_hndl.eq.0) then
          call set_gs_zavg_hndl(gs_avg_hndl,nelxy,ifld)
      endif

      nel = nelfld(ifld)
      n   = nx1*ny1*nz1*nel

      call copy(ua,bm1,n)              ! Set the averaging weights
      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weights over columns


      do i=1,n                          ! ua = (w_j*u_j)/( sum_i w_i)
         ua(i,1,1,1) = bm1(i,1,1,1)*u(i,1,1,1)/ua(i,1,1,1)
      enddo

      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weighted values


      return
      end
c-----------------------------------------------------------------------
      subroutine set_gs_zavg_hndl(gs_avg_hndl,nelxy,ifld)

c     Set the z-average handle

      include 'SIZE'
      include 'TOTAL'

      integer gs_avg_hndl,e,ex,ey,ez,eg

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /c_is1/ glo_num(lx1,ly1,lz1,lelv)
      integer*8 glo_num,ex_g


      nel = nelfld(ifld)
      do e=1,nel
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxy,1,1)

         ex_g = ex       ! Ensure int*8 promotion
         do k=1,nz1      ! Enumerate points in the x-y plane
            do j=1,ny1
            do i=1,nx1
               glo_num(i,j,k,e) = i+nx1*(j-1) + nx1*ny1*(ex_g-1)
            enddo
            enddo
         enddo

      enddo

      n = nel*nx1*ny1*nz1

      call gs_setup(gs_avg_hndl,glo_num,n,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine tke2_k10(fname1127,fname2127) ! simple average of files
c     computes tke 
      include 'SIZE'
      include 'TOTAL'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      real uk(lt),vk(lt),wk(lt)
      integer gs_avg_hndl
      save    gs_avg_hndl
      data    gs_avg_hndl / 0 /
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
     $    ,pa (lx1,ly2,lz2)     ,pb (lx1,ly1,lz2)
      real vxa(lx1*ly1*lz1*lelt)
      real vya(lx1*ly1*lz1*lelt)
      real vza(lx1*ly1*lz1*lelt)
      real pma(lx1*ly1*lz1*lelt)
      real tza(lx1*ly1*lz1*lelt)
      character*127 fname1127
      character*127 fname2127
      call auto_averager(fname1127)
      call copy(uk,vx) 
      call copy(vk,vy) 
      call copy(wk,vz) 
      call col2(uk,uk,lt)
      call col2(vk,vk,lt)
      call col2(wk,wk,lt)
      call auto_averager(fname2127)
      call sub2(vx,uk,lt)
      call sub2(vy,vk,lt)
      call sub2(vz,wk,lt)
      nelxy = 3068
      ifld  = 1
      call z_avg(vxa,vx,gs_avg_hndl,nelxy,ifld)
      call z_avg(vya,vy,gs_avg_hndl,nelxy,ifld)
      call z_avg(vza,vz,gs_avg_hndl,nelxy,ifld)
      call outpost(vxa,vya,vza,vxa,vya,'   ')
      if (nid.eq.0) write(6,*) 'Space averaged results written'
      call outfld2d_k10(vxa,vya,vza,pma,tza,nelxy,'tke',1)
      if (nid.eq.0) write(6,*) 'Time-space 2D results written'

      return
      end
c-----------------------------------------------------------------------
      subroutine heseig(Ac,n,m,wr,wi)
      include 'SIZE'
      include 'TOTAL'
c     Ac is matrix of size nxn and Bc is the matrix mxm from Ac
      integer i,ns,info,f,n,m,ind1,flag
      parameter (lbw=4*lx1*ly1*lz1*lelv)
      common /bigw/ bw(lbw)
      real Ac(1),B(m,m),wr(m),wi(m),z,H(1)

      call rzero(B,m**2)

      do i=1,m
      do j=1,m
         ind1 = (i-1)*n+j
         B(j,i) = Ac(ind1)
      enddo
      enddo

      call DHSEQR('E','N',m,1,m,B,m,wr,wi,z,1,bw,lbw,info)

      return
      end
c-----------------------------------------------------------------------
      subroutine hescon(Ac,n,m,rcond,normt)
      include 'SIZE'
      include 'TOTAL'
c     Ac is matrix of size nxn and Bc is the matrix mxm from Ac
      integer i,ns,info,f,n,m,ind1,flag
      real kwk1(4*m)
      integer kwk2(m)
      real Ac(1),B(m,m),piv(m),anorm,rcond
      character*1 normt

      call rzero(B,m**2)

      do i=1,m
      do j=1,m
         ind1 = (i-1)*n+j
         B(j,i) = Ac(ind1)
      enddo
      enddo

      anorm =  DLANGE( normt, m, m, B, m, kwk1 )
      call DGETRF( m, m, B, m, piv, info )
      call DGECON( normt, m, B, m, anorm, rcond,kwk1,kwk2,info)

      return
      end
c-----------------------------------------------------------------------
      subroutine getcondmax(Ac,n,m,normt,cond,maxeig,mineig)
      include 'SIZE'
      include 'TOTAL'
      real cond,eigr(m),eigi(m),j,k,n1,n2,n3,eigmag(m)
      integer i,m,n
      real invcond,maxeig,mineig
      character*1 normt
      real Ac(1)
        

      call heseig(Ac,n,m,eigr,eigi) !gets the eigenvalues
      call hescon(Ac,n,m,invcond,normt)  !gets the inverse of condition number
      maxeig = -1.e7
      mineig = 1.e7
      do i=1,m
         eigmag(i) = sqrt(eigr(i)**2+eigi(m)**2)
         if (eigmag(i).gt.maxeig) maxeig=eigmag(i)
         if (eigmag(i).lt.mineig) mineig=eigmag(i)
      enddo       

       cond = 0.
       if (invcond.ne.0) cond = 1./invcond


      return
      end
c-----------------------------------------------------------------------
      subroutine kk10fixcurs
      include 'SIZE'
      include 'TOTAL'
c     This routine fixes the curved edges post smoothing
      common /ctmp0/ zg(3)
      real w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelt)
      integer i,j,f,k,e,eg,nedge,n,nfaces
      integer ke(3,12)
      real xyz(3,3)
      real xd(lx1*ly1*lz1),yd(lx1*ly1*lz1)
      real zd(lx1*ly1*lz1)
      real nx,ny,nz,nxv,nyv
      save    ke
      data    ke /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1
     $           , 19,20,21,   21,24,27,   27,26,25,   25,22,19
     $           ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /
       integer kin(7)
       real tol

c     START BY LOOPING OVER EACH ELEMENT AND THEN OVER EACH EDGE
      nedge = 4 + 8*(ldim-2)
      nfaces = 2**ldim

      n      = nx1*ny1*nz1*nelv
      call rone  (v1mask,n)
      do e=1,nelv                ! fill mask where bc is periodic
      do f=1,nfaces              ! so we don't translate periodic bcs (z only)
       if (cbc(f,e,1).ne.'E  ') call facev (v1mask,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo
      call dsop(v1mask,'*  ',nx1,ny1,nz1)

      do e=1,nelv
       do j=1,nedge
        call rzero(xyz,9)
        do i=1,3
         xyz(1,i) = xm1(ke(i,j),1,1,e)
         xyz(2,i) = ym1(ke(i,j),1,1,e)
         if (ldim.eq.3) xyz(3,i) = zm1(ke(i,j),1,1,e)
        enddo

        if (v1mask(ke(2,j),1,1,e).gt.0.5) then
         call kk10fixedg2(xyz,nx,ny,nz)
         xm1(ke(2,j),1,1,e) = nx
         ym1(ke(2,j),1,1,e) = ny
         if (ldim.eq.3) zm1(ke(2,j),1,1,e) = nz
        endif
       enddo

c      fix face center and body center
        zg(1) = -1
        zg(2) = 0
        zg(3) = 1
        call copy(xd,xm1(1,1,1,e),lx1*ly1*lz1)
        call copy(yd,ym1(1,1,1,e),lx1*ly1*lz1)
        if (ldim.eq.3) call copy(zd,zm1(1,1,1,e),lx1*ly1*lz1)
        call gh_face_extend(xd,zg,3,2,w1,w2)
        call gh_face_extend(yd,zg,3,2,w1,w2)
        if (ldim.eq.3) call gh_face_extend(zd,zg,3,2,w1,w2)
        kin(1) = 5
        kin(2) = 11
        kin(3) = 13
        kin(4) = 15
        kin(5) = 17
        kin(6) = 23
        kin(7) = 14
        do i=1,(ldim-2)*6+1
         xm1(kin(i),1,1,e) = xd(kin(i))
         ym1(kin(i),1,1,e) = yd(kin(i))
         if (ldim.eq.3) zm1(kin(i),1,1,e) = zd(kin(i))
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine kk10fixedg2(xyz,nx,ny,nz)

      real xyz(3,0:2) ! Coordinates

      real x0(3),x1(3),x2(3),v2(3),v1(3),l02,l01,n1i,u1(3),u2(3)
      real nx,ny,nz
      real h,tol,h2

      call copy(x0,xyz(1,0),3)
      call copy(x1,xyz(1,1),3)
      call copy(x2,xyz(1,2),3)
                             !                  ^ x2
      call sub3(v1,x1,x0,3)  ! v1=x1-x0        /  \ v2     p1 = x0+v2*dot
      call sub3(v2,x2,x0,3)  ! v2=x2-x0    x1 .<-- x0
                             !                  v1
      l02 = vlsc2(v2,v2,3)   ! || x2-x0 ||
      if (l02.gt.0) then
         l02 = sqrt(l02)
         scl = 1./l02
         call cmult2(u2,v2,scl,3)  ! Unit tangent
      endif

      l01 = vlsc2(v1,v1,3)   ! || x1-x0 ||
      if (l01.gt.0) then
         l01 = sqrt(l01)
         scl = 1./l01
         call cmult2(u1,v1,scl,3)  ! Unit tangent
      endif

      dot = vlsc2(v1,u2,3)

      if (dot.le.0) then
         write(6,*) 'ERROR 1 IN SHIFT - YOU SHOULD ABORT'
         return
      elseif (dot.gt.l02) then
         write(6,*) 'ERROR 2 IN SHIFT - YOU SHOULD ABORT'
         return
      endif
         h = 0.
         tol = 1.e-8
      do i=1,3
         p1i = x0(i) + u2(i)*dot   ! Projection of x1 onto [x0,x2]
         n1i = x1(i) - p1i         ! Normal vector
         h = h + n1i**2
         h2 = h2+ (x2(i)-x0(i))**2
         xmi = 0.5*(x0(i)+x2(i))
         x1(i) = xmi + n1i         ! X1 point shifted to be centered at midpoint
      enddo
         if (h.le.h2*tol) then
            x1(1) = 0.5*(x0(1)+x2(1))
            x1(2) = 0.5*(x0(2)+x2(2))
            x1(3) = 0.5*(x0(3)+x2(3))
         endif


      nx = x1(1)
      ny = x1(2)
      nz = x1(3)

      return
      end
c-----------------------------------------------------------------------
      subroutine xx_avg(ua,u,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
      include 'SIZE'
      include 'TOTAL'

      real u (lx1,ly1,lz1,lelt)
      real ua(lx1,ly1,lz1,lelt)

      integer gs_avg_hndl,e,ex,ey,ez,eg

      if (gs_avg_hndl.eq.0) then
          call set_gs_xxavg_hndl(gs_avg_hndl,nelxx,nelyy*nelzz,ifld)
      endif

      nel = nelfld(ifld)
      n   = nx1*ny1*nz1*nel

      call copy(ua,bm1,n)              ! Set the averaging weights
      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weights over columns


      do i=1,n                          ! ua = (w_j*u_j)/( sum_i w_i)
         ua(i,1,1,1) = bm1(i,1,1,1)*u(i,1,1,1)/ua(i,1,1,1)
      enddo

      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weighted values


      return
      end
c-----------------------------------------------------------------------
      subroutine set_gs_xxavg_hndl(gs_avg_hndl,nelxx,nelyz,ifld)

c     Set the z-average handle

      include 'SIZE'
      include 'TOTAL'

      integer gs_avg_hndl,e,ex,ey,ez,eg,nelyz,nelxx

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /c_is1/ glo_num(lx1,ly1,lz1,lelv)
      integer*8 glo_num,ex_g


      nel = nelfld(ifld)
      do e=1,nel
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxx,nelyz,1)

         ex_g = ey       ! Ensure int*8 promotion
         do k=1,nz1      ! Enumerate points in the x-y plane
            do j=1,ny1
            do i=1,nx1
               glo_num(i,j,k,e) = j+ny1*(k-1) + ny1*nz1*(ex_g-1)
            enddo
            enddo
         enddo

      enddo

      n = nel*nx1*ny1*nz1

      call gs_setup(gs_avg_hndl,glo_num,n,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
