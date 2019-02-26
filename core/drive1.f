c-----------------------------------------------------------------------
      subroutine nek_init(comm_out)
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'

C     used scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs automatically
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

      integer comm_out
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
  
      common /rdump/ ntdump

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 

      common /c_is1/ glo_num(lx1 * ly1 * lz1, lelt)
      common /ivrtx/ vertex((2 ** ldim) * lelt)
      integer*8 glo_num, ngv
      integer vertex

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8) 
      ! set word size for LOGICAL
      lsize = sizeof(ltest) 
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm()
      nekcomm  = intracomm
      comm_out = nekcomm
      call iniproc()

      etimes = dnekclock()
      istep  = 0

      call opcount(1)

      call initdim         ! Initialize / set default values.
      call initdat
      call files

      etime = dnekclock()
      call readat          ! Read .rea +map file

      etims0 = dnekclock_sync()
      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
 12      format(1X,A,4I12,/,/)
      endif 

      call setvar          ! Initialize most variables

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call fix_geom
      call geom_reset(1)    ! recompute Jacobians, etc.

      call vrdsmsh          ! verify mesh topology
      call mesh_metrics     ! print some metrics

      call setlog  ! Initalize logical flags

      if (ifneknekc) call neknek_setup
      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0 .or. nsteps.ne.0) 
     $   call geneig(igeom) ! eigvals for tolerances

      call dg_setup ! Setup DG, if dg flag is set.

      if (ifflow.and.iftran) then ! Init pressure solver 
         if (fintim.ne.0 .or. nsteps.ne.0) call prinit
      endif

      if(ifcvode) call cv_setsize

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

#ifdef CMTNEK
        call nek_cmt_init
#endif

      call setics
      call setprop

      if (instep.ne.0) then
         if (ifneknekc) call neknek_exchange
         if (ifneknekc) call chk_outflow

         if (nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif

      call setprop      ! call again because input has changed in userchk

      if (ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0            ! Set perturbation field count to 0 for baseline flow

      call in_situ_init()

      call time00       !     Initalize timers to ZERO
      call opcount(2)

      ntdump=0
      if (timeio.ne.0.0) ntdump = int( time/timeio )

      tinit = dnekclock_sync() - etimes
      if (nio.eq.0) then
        write (6,*) ' '
        if (time.ne.0.0) write (6,'(a,e14.7)') ' Initial time:',time
        write (6,'(a,g13.5,a)') 
     &     ' Initialization successfully completed ', tinit, ' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'NEKNEK'
      logical if_ms_multdt
      common /int_logical_ms/ if_ms_multdt

      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      if(ifsync) isyc=1
      itime = 0
#ifdef TIMER
      itime = 1
#endif
      call nek_comm_settings(isyc,itime)

      ! start measurements
      call nek_comm_startstat()
      dtmp = dnekgflops()

      istep  = 0
      msteps = 1

      do kstep=1,nsteps,nss_ms
         if (if_ms_multdt) then 
           call nek__multi_advance_dt(kstep,nss_ms)
         else
           call nek__multi_advance(kstep,msteps)
         endif
         if(kstep.ge.nsteps) lastep = 1
         call check_ioinfo  
         call set_outfld
         etime1 = dnekclock()
         call userchk
         tuchk = tuchk + dnekclock()-etime1
         call prepost (ifoutfld,'his')
         call in_situ_check()
         if (mod(kstep,100).eq.0 ..and. lastep.eq.0) call runstat
         if (lastep .eq. 1) goto 1001
      enddo
 1001 lastep=1


      call nek_comm_settings(isyc,0)

      call comment

c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nio.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif


      RETURN
      END

c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      ntot = lx1*ly1*lz1*nelv

      call nekgsync

      call setup_convect(2) ! Save conv vel

      if (iftran) call settime
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

#ifdef CMTNEK
      if (nio.eq.0.and.istep.le.1) write(6,*) 'CMT branch active'
      call cmt_nek_advance
      return
#endif

      if (ifsplit) then   ! PN/PN formulation

         do igeom=1,ngeom

         if (ifneknekc .and. igeom.gt.2) then
            if (ifneknekm.and.igeom.eq.3) call neknek_setup
            call neknek_exchange
         endif

         ! call here before we overwrite wx 
         if (ifheat .and. ifcvode) call heat_cvode (igeom)   

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif

         if (ifheat) call heat (igeom)

         if (igeom.eq.2) then  
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif

         if (ifflow)          call fluid    (igeom)
         if (ifmvbd)          call meshv    (igeom)
         if (igeom.eq.ngeom.and.filterType.eq.1)
     $                        call q_filter(param(103))

         enddo

      else                ! PN-2/PN-2 formulation
         call setprop
         do igeom=1,ngeom

            if (ifneknekc .and. igeom.gt.2) then
              if (ifneknekm.and.igeom.eq.3) call neknek_setup
              call neknek_exchange
            endif

            ! call here before we overwrite wx 
            if (ifheat .and. ifcvode) call heat_cvode (igeom)   

            if (ifgeom) then
               if (.not.ifrich) call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif
            if (igeom.eq.ngeom.and.filterType.eq.1)
     $         call q_filter(param(103))
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TOTAL'

      if(instep.ne.0) call runstat

c      if (ifstrs) then
c         call fgslib_crs_free(xxth_strs) 
c      else
c         call fgslib_crs_free(xxth(1))
c      endif

      call in_situ_end()
      call exitt0()

      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance(kstep,msteps)

      include 'SIZE'
      include 'TOTAL'

      do i=1,msteps
         istep = istep+i
         call nek_advance

         if (ifneknekc) then 
            call neknek_exchange
            call bcopy
            call chk_outflow
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance_dt(kstep,msteps)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      real prsav(lx2*ly2*lz2*lelv,0:10)
      real prlagdt(lx2*ly2*lz2*lelv)

      real vxyzsav(lx1*ly1*lz1*lelv,ldim,0:10)
      real vxyzd(lx1*ly1*lz1*lelv,ldim,0:1)

      real tsav(lx1*ly1*lz1*lelt,ldimt,0:10)
      real td(lx1*ly1*lz1*lelt,ldimt,0:1)

      real timsav

      ntotv = lx1*ly1*lz1*nelv
      ntotp = lx2*ly2*lz2*nelv
      ntott = lx1*ly1*lz1*nelt
!     Copy u,v,w,pr,t at t^{n-1}
!     Copy lagging arrays for BDF terms 
      if (ifflow) then
        call copy(prsav(1,0),pr,ntotp)        !p solution
        call copy(vxyzsav(1,1,0),vx,ntotv)    !vx solution
        call copy(vxyzsav(1,2,0),vy,ntotv)    !vy solution
        if (ldim.eq.3) call copy(vxyzsav(1,ldim,0),vz,ntotv)
        do j=1,2
          call copy(vxlagdt(1,1,1,1,j),vxlag(1,1,1,1,j),ntotv) !vxlag
          call copy(vylagdt(1,1,1,1,j),vylag(1,1,1,1,j),ntotv) !vylag
         if (ldim.eq.3)
     $    call copy(vzlagdt(1,1,1,1,j),vzlag(1,1,1,1,j),ntotv) !vzlag
        enddo
        if (.not.ifsplit) call copy(prlagdt,prlag,ntotp)       !prlag

        call copy(abx1dt,abx1,ntotv)                           !nl term
        call copy(aby1dt,aby1,ntotv)
        if (ldim.eq.3) call copy(abz1dt,abz1,ntotv)
        call copy(abx2dt,abx2,ntotv)
        call copy(aby2dt,aby2,ntotv)
        if (ldim.eq.3) call copy(abz2dt,abz2,ntotv)
      endif
      if (ifheat) then
        call copy(tsav(1,1,0),t(1,1,1,1,1),ntott) !t solution
        do j=1,2
          call copy(tlagdt(1,1,1,1,j,1),tlag(1,1,1,1,j,1),ntott) !t lag
        enddo
        call copy(vgradt1dt(1,1,1,1,1),vgradt1(1,1,1,1,1),ntott) !u.del t
        call copy(vgradt2dt(1,1,1,1,1),vgradt2(1,1,1,1,1),ntott)
      endif

      if (istep.eq.0.and.ifneknekc) then
        call neknek_exchange
        call bcopy_only
      endif

      !exchange information at t^{n-1} and save it.. we need this for 
      !boundary condition during corrector iterations
      if (ifflow) then
        do j=1,ldim
         call copy(vxyzd(1,j,0),valint(1,1,1,1,j),ntotv)
        enddo
      endif
      if (ifheat) then
        call copy(td(1,1,0),valint(1,1,1,1,ldim+2),ntott)
      endif
      
      iorigstep = istep
      timsav = time

      do i=1,msteps
        iss_ms = i !multisession sub-step
        istep = istep+1
c       calculate appropriate extrapolation coefficients
        rcoeff = i*1./msteps
        if (msteps.eq.1) rcoeff = itstepratio*1.
        call bdr_data_extr(rcoeff)
c       solve with ngeom=2 
        call nek_advance_ms(1,2,1)
c       save u,v,w,p,t at each sub-step
        if (ifflow) then
         call copy(vxyzsav(1,1,i),vx,ntotv)
         call copy(vxyzsav(1,2,i),vy,ntotv)
         if (ldim.eq.3) call copy(vxyzsav(1,ldim,i),vz,ntotv)
         call copy(prsav(1,i),pr,ntotp)
        endif
        if (ifheat) then
         call copy(tsav(1,1,i),t(1,1,1,1,1),ntott)
        endif
      enddo

      call neknekgsync()

c     transfer latest solutions and save them. These will also 
c     be used during corrector iterations (updated after each iteration)
      if (ifflow) then
        call neknek_xfer_fld(vx,vxyzd(1,1,1))
        call neknek_xfer_fld(vy,vxyzd(1,2,1))
        if (ldim.eq.3) call neknek_xfer_fld(vz,vxyzd(1,ldim,1))
      endif
      if (ifheat) then
        call neknek_xfer_fld(t(1,1,1,1,1),td(1,1,1))
      endif

cc    Schwarz iterations
      do igeom=3,ngeom
c       Restor u,v,w,pr,t at t^{n-1} along with lagging arrays
        if (ifflow) then
          call copy(vx,vxyzsav(1,1,0),ntotv)
          call copy(vy,vxyzsav(1,2,0),ntotv)
          if (ldim.eq.3) call copy(vz,vxyzsav(1,ldim,0),ntotv)
          call copy(pr,prsav(1,0),ntotp)
          do j=1,2
            call copy(vxlag(1,1,1,1,j),vxlagdt(1,1,1,1,j),ntotv)
            call copy(vylag(1,1,1,1,j),vylagdt(1,1,1,1,j),ntotv)
            call copy(vzlag(1,1,1,1,j),vzlagdt(1,1,1,1,j),ntotv)
          enddo
          if (.not.ifsplit) call copy(prlag,prlagdt,ntotp)
          call copy(abx1,abx1dt,ntotv)
          call copy(aby1,aby1dt,ntotv)
          if (ldim.eq.3) call copy(abz1,abz1dt,ntotv)
          call copy(abx2,abx2dt,ntotv)
          call copy(aby2,aby2dt,ntotv)
          if (ldim.eq.3) call copy(abz2,abz2dt,ntotv)
        endif !ifflow
        if (ifheat) then
          call copy(t(1,1,1,1,1),tsav(1,1,0),ntott)
          call copy(vgradt1(1,1,1,1,1),vgradt1dt(1,1,1,1,1),ntott)
          call copy(vgradt2(1,1,1,1,1),vgradt2dt(1,1,1,1,1),ntott)
          do j=1,2
            call copy(tlag(1,1,1,1,j,1),tlagdt(1,1,1,1,j,1),ntott)
          enddo
        endif !ifheat

        istep = iorigstep
        time = timsav
        do i=1,msteps
          istep = istep+1
          c1 = i*1./(msteps*1.)
          c0 = 1.-c1
          if (ifflow) then
            do j=1,ldim
              call add3s2(valint(1,1,1,1,j),vxyzd(1,j,0),vxyzd(1,j,1),
     $                    c0,c1,ntotv)
            enddo
          endif
          if (ifheat) then
             call add3s2(valint(1,1,1,1,ldim+2),td(1,1,0),td(1,1,1),
     $                   c0,c1,ntott)
          endif

          igeomo = igeom
          call nek_advance_ms(1,igeomo,igeom-1)

          if (ifflow) then
            call copy(vxyzsav(1,1,i),vx,ntotv)
            call copy(vxyzsav(1,2,i),vy,ntotv)
            if (ldim.eq.3) call copy(vxyzsav(1,ldim,i),vz,ntotv)
            call copy(prsav(1,i),pr,ntotp)
          endif !ifflow
          if (ifheat) then
            call copy(tsav(1,1,i),t(1,1,1,1,1),ntott)
          endif !ifheat
        enddo !msteps i.e. substeps
        call neknekgsync()
 
        if (ifflow) then
          call neknek_xfer_fld(vx,vxyzd(1,1,1))
          call neknek_xfer_fld(vy,vxyzd(1,2,1))
          if (ldim.eq.3) call neknek_xfer_fld(vz,vxyzd(1,ldim,1))
        endif
        if (ifheat) then
          call neknek_xfer_fld(t(1,1,1,1,1),td(1,1,1))
        endif
      enddo !igeom loop

      call neknek_exchange
      call bcopy_only

      if (itstepratio.gt.1) then
        n=0
        do j=itstepratio-1,itstepratio-2,-1
           n=n+1
           if (ifflow) then
              do k=1,ldim
               call neknek_xfer_fld(vxyzsav(1,k,j),vxyzd(1,k,1))
               if (msteps.eq.1) call copy(bdrylg(1,k,n),
     $                                   vxyzd(1,k,1),ntotv)
              enddo 
           endif !ifflow
           if (ifheat) then
                 call neknek_xfer_fld(tsav(1,1,j),td(1,1,1))
                 if (msteps.eq.1) call copy(bdrylg(1,ldim+2,n),
     $                                   td(1,1,1),ntott)
           endif !ifheat
        enddo !j=itstepratio-1,itstepratio-2,-1
      endif !itstepratio.gt.1

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_advance_ms(igeomstart,igeomend,igeomskip)
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      ntot = lx1*ly1*lz1*nelv

      call nekgsync

      if (igeomstart.ne.1) goto 1002
      call setup_convect(2) ! Save conv vel

      if (iftran) call settime
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

 1002 continue

      if (ifsplit) then   ! PN/PN formulation

         do igeom=igeomstart,igeomend,igeomskip

c         if (ifneknekc .and. igeom.gt.2) then
c            if (ifneknekm.and.igeom.eq.3) call neknek_setup
c            call neknek_exchange
c         endif

         ! call here before we overwrite wx 
         if (ifheat .and. ifcvode) call heat_cvode (igeom)

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif

         if (ifheat) call heat (igeom)

         if (igeom.ge.2) then
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif

         if (ifflow)          call fluid    (igeom)
         if (ifmvbd)          call meshv    (igeom)
         if (igeom.eq.ngeom.and.filterType.eq.1)
     $                        call q_filter(param(103))

         enddo

      else                ! PN-2/PN-2 formulation
c        call exitti('Pn-Pn-2 currently disabled$',lelt)
         call setprop
         do igeom=igeomstart,igeomend,igeomskip

c           if (ifneknekc .and. igeom.gt.2) then
c             if (ifneknekm.and.igeom.eq.3) call neknek_setup
c             call neknek_exchange
c           endif

            ! call here before we overwrite wx 
            if (ifheat .and. ifcvode) call heat_cvode (igeom)

            if (ifgeom) then
               if (.not.ifrich) call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif
            if (igeom.eq.ngeom.and.filterType.eq.1)
     $         call q_filter(param(103))
         enddo
      endif
      igeom = igeomend

      return
      end
c-----------------------------------------------------------------------
