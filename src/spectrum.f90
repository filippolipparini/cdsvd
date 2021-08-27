program cd_spectrum
implicit none
!
! given a set of rotational strenght tensors and the standard deviation of
! a gaussian shape function, build a CD spectrum of the molecule in a random
! orientation.
!
! This is done by generating a random rotation matrix and using it to rotate
! the rotational strength tensors; then, the zz component is used as an 
! intensity.
!
!integer, allocatable ::
real*8,  allocatable :: rs(:,:,:), rsrot(:,:,:), spectrum(:,:), rszz(:), freq(:), fave(:)
real*8,  allocatable :: sp_ave(:)
!
integer  :: nstates, nspectra, npoint, lmin, lmax, it, i, j, istate, ispectrum, idata
integer  :: ijunk, lfname
real*8   :: kmat(3,3), umat(3,3), tmp(3,3), fwhm, xx, yy, zz, rand, norm
real*8   :: xrand, yrand, zrand, thrand, pi
real*8   :: eev, fac, sigma, gaussian
logical  :: done
real*8, parameter :: zero=0.0d0, pt5=0.50d0, one=1.0d0, two=2.0d0, four=4.0d0, tol=1.0d-12
character*20 cform, args(5), fname
!
! read in some information:
!
if (iargc().ne.5) then
  write(6,*) ' please run this program as follows: '
  write(6,*) ' ./spectrum lmin lmax nspectra, fwhm, idata'
  write(6,*) '   lmin/max: minimum and maximum wavelengths'
  write(6,*) '   nspectra: how many spectra you want to generate'
  write(6,*) '   fwhm:     fwhm for gaussian convolution'
  write(6,*) '   idata:    1 for full tensor, 2 for mu.m, 3 for mu.q'
  stop
end if
do i = 1, 5
  call getarg(i,args(i))
end do
read(args(1),*) lmin
read(args(2),*) lmax
read(args(3),*) nspectra
read(args(4),*) fwhm
read(args(5),*) idata
npoint = lmax - lmin + 1
if (idata .eq. 1) then
  fname = 'tensors.dat'
else if (idata .eq. 2) then
  fname = 'tensor-Rm.dat'
else if (idata .eq. 3) then
  fname = 'tensor-Rq.dat'
else
  write(6,*) ' invalid selection'
  stop
end if
!
fac   = sqrt(two*log(two))
sigma = fwhm/fac
pi    = four*atan(one)
!
allocate (spectrum(npoint,nspectra))
!
open (unit=10,file=fname,form='formatted',access='sequential')
read(10,*) nstates
allocate (freq(nstates), rs(3,3,nstates), rsrot(3,3,nstates), rszz(nstates))
do i = 1, nstates
  read(10,*) freq(i)  ! excitation energy in eV
end do
do i = 1, nstates
  do j = 1, 3
    read(10,*) ijunk, xx, yy, zz
    rs(j,1,i) = xx
    rs(j,2,i) = yy
    rs(j,3,i) = zz
  end do
end do
close (10)
!
spectrum = zero
!
do ispectrum = 1, nspectra
  kmat = zero
!
! generate a random rotation axis:
!
  call random_number(xrand)
  call random_number(yrand)
  call random_number(zrand)
  fac   = sqrt(xrand**2 + yrand**2 + zrand**2)
  xrand = xrand/fac
  yrand = yrand/fac
  zrand = zrand/fac
!
! generate a random rotation angle:
!
  call random_number(thrand)
  thrand = thrand*two*pi
!
! form the skew-symmetric matrix kappa
!
  kmat(1,2) = xrand*thrand
  kmat(1,3) = yrand*thrand
  kmat(2,3) = zrand*thrand
  kmat(2,1) = - kmat(1,2)
  kmat(3,1) = - kmat(1,3)
  kmat(3,2) = - kmat(2,3)
  umat = zero
  umat(1,1) = one
  umat(2,2) = one
  umat(3,3) = one
!
  done = .false.
!
! generate the unitary transformation
!
  tmp = kmat
  it  = 1
  do while (.not. done)
    norm = dot_product(tmp(:,1),tmp(:,1)) + dot_product(tmp(:,2),tmp(:,2)) + dot_product(tmp(:,3),tmp(:,3))
    norm = sqrt(norm)
    done = norm.lt.tol
    it   = it + 1
    fac  = one/float(it)
    umat = umat + tmp
    tmp  = fac*matmul(tmp,kmat)
  end do
!
! rotate all the states:
!
  do istate = 1, nstates
    tmp   = matmul(transpose(umat),rs(:,:,istate))
    rsrot(:,:,istate) = matmul(tmp,umat)
    rszz(istate) = rsrot(3,3,istate)
  end do
!
! now assemble the spectrum by convoluting the signals with gaussian line shapes:
!
  do i = 1, npoint
    eev = 1240.d0/float(lmin+i-1)
    do istate = 1, nstates
      spectrum(i,ispectrum) = spectrum(i,ispectrum) + rszz(istate)*gaussian(eev,freq(istate),sigma)
    end do
  end do
!
end do
!
! write the spectrum on a file:
!
open (unit=11,file='data.txt',form='formatted',access='sequential')
write(cform,'("(I6,",I4,"F16.8)")') nspectra
do i = 1, npoint
  write(11,cform) lmin + i - 1, spectrum(i,:)
end do
close (11)
!
! compute the isotropic spectrum:
!
allocate (fave(nstates),sp_ave(npoint))
do istate = 1, nstates
  fave(istate) = (rs(1,1,istate) + rs(2,2,istate) + rs(3,3,istate))/3.0d0
end do
!
sp_ave = zero
do i = 1, npoint
  eev = 1240.d0/float(lmin+i-1)
  do istate = 1, nstates
    sp_ave(i) = sp_ave(i) + fave(istate)*gaussian(eev,freq(istate),sigma)
  end do
end do
!
open (unit=12,file='iso.dat',form='formatted',access='sequential')
do i = 1, npoint
  write(12,*) lmin + i - 1, sp_ave(i)
end do
close(12)

end program cd_spectrum
!
! function gaussian(x,center,fwhm)
!
real*8 function gaussian(x,center,sigma)
implicit none
real*8            :: x, center, sigma
!
real*8            :: pi, fnorm, fexp
real*8, parameter :: one=1.0d0, two=2.0d0, four=4.0d0
!
pi    = four*atan(one)
fnorm = one/(sigma*sqrt(two*pi))
!
fexp  = - (x - center)**2
fexp  = fexp / (two*sigma*sigma)
gaussian = fnorm*exp(fexp)
return
end function gaussian
