program cdsvd
  implicit none
!
  integer :: i, j, nrow, ncol, lwork, info, ios, nprint, nsing, k, klim
  real*8  :: xx(1), lw(1), error, x, y, z
  real*8, allocatable :: freq(:), a(:,:), u(:,:), vt(:,:), s(:), work(:), cds(:,:), cdave(:,:), oldave(:)
  character*20 cform, args(2)
!
  if (iargc().ne.2) then
    write(6,*) ' please run this program as follows: '
    write(6,*) ' ./cdsvd nrow ncol'
    stop
  end if
!
! gather information on the size of the matrix:
!
  do i = 1, 2
    call getarg(i,args(i))
  end do
  read(args(1),*) nrow
  read(args(2),*) ncol
!
! read the matrix
!
  allocate (a(nrow,ncol),freq(nrow))
  open (unit=1,file='data.txt',form='formatted',access='sequential',iostat=ios)
  if (ios.ne.0) then
    write(6,*) ' error opening the data file.'
    write(6,*) ' the file must be named "data.txt"'
    stop
  end if

  do i = 1, nrow
    read(1,*) freq(i), a(i,1:ncol)
  end do
  close(1)
!
  allocate (oldave(nrow))
  oldave = 0.0d0
  do j = 1, ncol
    oldave = oldave + a(:,j)
  end do
  oldave = oldave/float(ncol)
!
! set up the svd. get memory...
!
  call dgesvd('a','a',nrow,ncol,xx,nrow,xx,xx,nrow,xx,ncol,lw,-1,info)
  lwork = int(lw(1))
!
  nsing = min(nrow,ncol)
  allocate (s(min(nrow,ncol)), u(nrow,nrow), vt(ncol,ncol), work(lwork))
!
  call dgesvd('a','a',nrow,ncol,a,nrow,s,u,nrow,vt,ncol,work,lwork,info)
!
  write(6,*) ' singular values of the matrix:', nsing
  write(6,'(10f14.8)') s
!
! assemble the singular cd spectra:
!
  allocate (cds(nrow,ncol))
  do j = 1, ncol
    do k = 1, nsing
      do i = 1, nrow
        cds(i,k) = cds(i,k) + u(i,k)*vt(k,j)
      end do
    end do
  end do
!
  allocate (cdave(nrow,nsing))
  cdave = 0.0d0
  do klim = 1, nsing
    do k = 1, klim
      cdave(:,klim) = cdave(:,klim) + s(k)*cds(:,k)
    end do
    cdave(:,klim) = cdave(:,klim)/float(klim)
  end do
!
! safety check:
!
  oldave = oldave - cdave(:,nsing)
  error = sqrt(dot_product(oldave,oldave))
  if (error.gt.1.0d-6) then
    write(6,*) ' Waring: SVD error is ', error
  end if
  oldave = oldave + cdave(:,nsing)
!
  open (file='sing_val.dat',unit=10,form='formatted',access='sequential')
  do i = 1, nsing
    write(10,'(f16.8)') s(i)
  end do
  close (10)
  write(cform,'("(F16.8,",I4,"F16.8)")') nrow
  open (file='sing_cd.dat',unit=10,form='formatted',access='sequential')
  do i = 1, nrow
    write(10,cform) freq(i), cds(i,:)
  end do
  close (10)
  open (file='average.dat',unit=10,form='formatted',access='sequential')
  do i = 1, nrow
    write(10,'(2F16.8)') freq(i), oldave(i)
  end do
  close (10)
  deallocate(a,s,u,vt,work,freq,cdave,oldave,cds)
end
