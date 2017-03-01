c-----------------------------------------------------------------------
c NEKNEK
      subroutine get_valuesnbk10(ranvar,ranrec)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'mpif.h'

      integer nsize,i,j,k,n,id
      real ranvar
      real ranrec
      integer status(mpi_status_size)
C     Send interpolation values to the corresponding processors 
C     of remote session

      call neknekgsync()

      do id=0,npsend-1
         len=wdsize
         call mpi_irecv(ranrec,len,mpi_byte,id,id,intercomm,msg,ierr)
         call mpi_send(ranvar,len,mpi_byte,id,nid,intercomm,ierr)
         call mpi_wait (msg,status,ierr)
      enddo

      call neknekgsync()

      return
      end
C--------------------------------------------------------------------------
      subroutine get_valuesk10(ranvar,ranrec)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'mpif.h'

      integer nsize,i,j,k,n
      real ranvar
      real ranrec
      integer zzid,zzmat1(2)
      integer status(mpi_status_size)
      integer requ
c k10
      zzmat1(1) = 0
      zzmat1(2) = 1
c k10

      nsize = 1



      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

C     Send interpolation values to the corresponding processors 
C     of remote session

      call neknekgsync()
      do zzid=1,2

      if (idsess.eq.zzmat1(zzid)) then
      do id=0,npsend-1
         len=nsize*wdsize
         call mpi_send(ranvar,len,mpi_byte,id,nid,intercomm,ierr)
      enddo
      else
      do id=0,nprecv-1
       nrecv = nsize
       len=nrecv*wdsize
        call mpi_recv (ranrec,len,mpi_byte,id,id,intercomm,status,ierr)
      enddo

      endif
      enddo
      call neknekgsync()

      return
      end
C--------------------------------------------------------------------------
