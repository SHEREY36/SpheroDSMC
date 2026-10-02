!> Histograms of the reduced speed c* = |c - U|/sqrt(2 T_tr/m) and the reduced
!> rotational energy x = E_rot/T_rot, accumulated over sample times.
module distributions
  use precision, only: dp, i8, pi
  use particles, only: gas_t
  implicit none
  private
  public :: histogram_t, init_histograms, sample_histograms, write_histograms

  type :: histogram_t
    integer :: nbins = 0
    real(dp) :: width = 0
    integer(i8), allocatable :: count(:)
    integer(i8) :: total = 0, overflow = 0
  end type histogram_t

contains

  subroutine init_histograms(speed, rot)
    type(histogram_t), intent(out) :: speed, rot
    call init_one(speed, 400, 0.02_dp)      ! c* in [0, 8)
    call init_one(rot, 800, 0.05_dp)        ! x in [0, 40)
  end subroutine init_histograms

  subroutine init_one(h, nbins, width)
    type(histogram_t), intent(out) :: h
    integer, intent(in) :: nbins
    real(dp), intent(in) :: width
    h%nbins = nbins
    h%width = width
    allocate(h%count(nbins))
    h%count = 0
  end subroutine init_one

  subroutine sample_histograms(g, ttr, trot, speed, rot)
    type(gas_t), intent(in) :: g
    real(dp), intent(in) :: ttr, trot
    type(histogram_t), intent(inout) :: speed, rot
    real(dp) :: mean(3), vth
    integer :: i
    mean = sum(g%c, dim=2) / g%n
    vth = sqrt(2 * ttr / g%mass)
    do i = 1, g%n
      call add(speed, norm2(g%c(:, i) - mean) / vth)
      call add(rot, g%erot(i) / trot)
    end do
  end subroutine sample_histograms

  subroutine add(h, x)
    type(histogram_t), intent(inout) :: h
    real(dp), intent(in) :: x
    integer :: b
    h%total = h%total + 1
    b = int(x / h%width) + 1
    if (b > h%nbins) then
      h%overflow = h%overflow + 1
    else
      h%count(b) = h%count(b) + 1
    end if
  end subroutine add

  !> Normalised densities with their Maxwellian references:
  !> P_M(c*) = 4 c*^2 exp(-c*^2)/sqrt(pi) and P_M(x) = exp(-x).
  subroutine write_histograms(prefix, speed, rot)
    character(len=*), intent(in) :: prefix
    type(histogram_t), intent(in) :: speed, rot
    integer :: u, b
    real(dp) :: x
    open(newunit=u, file=prefix//'.speed', status='replace')
    write(u, '(a,i0,a,i0)') '# reduced speed c* = |c-U|/sqrt(2T_tr/m); samples ', speed%total, &
                            '; beyond range ', speed%overflow
    write(u, '(a)') '# c*  P(c*)  P_M(c*)  count'
    do b = 1, speed%nbins
      x = (b - 0.5_dp) * speed%width
      write(u, '(3es16.8,i14)') x, speed%count(b) / (real(max(speed%total, 1_i8), dp) * speed%width), &
                                4 * x * x * exp(-x * x) / sqrt(pi), speed%count(b)
    end do
    close(u)
    open(newunit=u, file=prefix//'.rot', status='replace')
    write(u, '(a,i0,a,i0)') '# reduced rotational energy x = E_rot/T_rot; samples ', rot%total, &
                            '; beyond range ', rot%overflow
    write(u, '(a)') '# x  P(x)  P_M(x)  count'
    do b = 1, rot%nbins
      x = (b - 0.5_dp) * rot%width
      write(u, '(3es16.8,i14)') x, rot%count(b) / (real(max(rot%total, 1_i8), dp) * rot%width), &
                                exp(-x), rot%count(b)
    end do
    close(u)
  end subroutine write_histograms

end module distributions
