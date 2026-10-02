!> Minimal reader for uncompressed NumPy archives written by np.savez.
!>
!> An .npz file is a ZIP archive of "stored" (uncompressed) .npy members.  Each
!> member is a short text header followed by the raw array.  We support the
!> layout produced by tools/export_model.py: little-endian float64 ('<f8') or
!> int64 ('<i8') arrays in C order, with up to four dimensions.
!>
!> A C array of shape (n0, n1) is returned as a Fortran array a(n1, n0), so
!> each C row is a contiguous Fortran column.
module npz_io
  use precision, only: dp, i8
  implicit none
  private
  public :: npz_read, npz_real1, npz_real2, npz_int1, npz_scalar

  type :: member_t
    character(len=64) :: name = ''
    integer :: ndim = 0
    integer(i8) :: shape(4) = 1
    logical :: is_int = .false.
    real(dp), allocatable :: r(:)
    integer(i8), allocatable :: k(:)
  end type member_t

  type, public :: npz_t
    character(len=:), allocatable :: path
    integer :: n = 0
    type(member_t), allocatable :: m(:)
  end type npz_t

contains

  subroutine npz_read(path, z)
    character(len=*), intent(in) :: path
    type(npz_t), intent(out) :: z
    character(len=:), allocatable :: buf, header
    integer :: unit, ios, nlen, elen, method, flags
    integer(i8) :: fsize, pos, ds, hstart, hlen, nbytes, first
    type(member_t), allocatable :: grown(:)

    open(newunit=unit, file=path, access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('cannot open '//trim(path))
    inquire(unit=unit, size=fsize)
    allocate(character(len=fsize) :: buf)
    read(unit) buf
    close(unit)

    z%path = path
    allocate(z%m(32))
    pos = 1
    do while (pos + 30 <= fsize)
      if (u32(buf, pos) /= 67324752_i8) exit            ! 'PK\3\4': local file header
      flags = int(u16(buf, pos + 6))
      method = int(u16(buf, pos + 8))
      nlen = int(u16(buf, pos + 26))
      elen = int(u16(buf, pos + 28))
      if (method /= 0) call fail(trim(path)//' is compressed: write it with np.savez')
      if (btest(flags, 3)) call fail(trim(path)//' uses ZIP data descriptors')

      if (z%n == size(z%m)) then
        allocate(grown(2 * size(z%m)))
        grown(1:z%n) = z%m(1:z%n)
        call move_alloc(grown, z%m)
      end if
      z%n = z%n + 1
      associate (mem => z%m(z%n))
        mem%name = buf(pos + 30:pos + 29 + nlen)
        if (index(mem%name, '.npy') > 0) mem%name = mem%name(1:index(mem%name, '.npy') - 1)

        ds = pos + 30 + nlen + elen                        ! start of the .npy member
        if (buf(ds:ds + 5) /= char(147)//'NUMPY') call fail(trim(mem%name)//' is not .npy')
        if (ichar(buf(ds + 6:ds + 6)) == 1) then
          hlen = u16(buf, ds + 8)
          hstart = ds + 10
        else
          hlen = u32(buf, ds + 8)
          hstart = ds + 12
        end if
        header = buf(hstart:hstart + hlen - 1)
        call parse_header(header, mem)

        nbytes = 8_i8 * product(mem%shape(1:max(mem%ndim, 1)))
        if (mem%ndim == 0) nbytes = 8_i8
        first = hstart + hlen
        if (mem%is_int) then
          mem%k = transfer(buf(first:first + nbytes - 1), 0_i8, int(nbytes / 8))
        else
          mem%r = transfer(buf(first:first + nbytes - 1), 0.0_dp, int(nbytes / 8))
        end if
        pos = first + nbytes
      end associate
    end do
    if (z%n == 0) call fail(trim(path)//' holds no .npy members')
  end subroutine npz_read

  !> Read descr, fortran_order and shape from a .npy header dictionary.
  subroutine parse_header(h, mem)
    character(len=*), intent(in) :: h
    type(member_t), intent(inout) :: mem
    integer :: i, j, k
    character(len=:), allocatable :: dims

    i = index(h, "'descr':")
    j = index(h(i + 8:), "'") + i + 8
    select case (h(j:j + 2))
    case ('<f8')
      mem%is_int = .false.
    case ('<i8')
      mem%is_int = .true.
    case default
      call fail(trim(mem%name)//': only <f8 and <i8 arrays are supported, got '//h(j:j + 2))
    end select
    if (index(h, "'fortran_order': True") > 0) call fail(trim(mem%name)//' is Fortran-ordered')

    i = index(h, "'shape':")
    j = index(h(i:), '(') + i
    k = index(h(j:), ')') + j - 2
    dims = h(j:k)
    mem%ndim = 0
    mem%shape = 1
    do while (len_trim(dims) > 0)
      i = index(dims, ',')
      if (i == 0) i = len(dims) + 1
      if (len_trim(dims(1:i - 1)) > 0) then
        mem%ndim = mem%ndim + 1
        if (mem%ndim > 4) call fail(trim(mem%name)//' has more than four dimensions')
        read(dims(1:i - 1), *) mem%shape(mem%ndim)
      end if
      if (i >= len(dims)) exit
      dims = dims(i + 1:)
    end do
  end subroutine parse_header

  integer function find(z, name) result(idx)
    type(npz_t), intent(in) :: z
    character(len=*), intent(in) :: name
    do idx = 1, z%n
      if (trim(z%m(idx)%name) == name) return
    end do
    call fail('array "'//name//'" not found in '//z%path)
  end function find

  subroutine npz_real1(z, name, a)
    type(npz_t), intent(in) :: z
    character(len=*), intent(in) :: name
    real(dp), allocatable, intent(out) :: a(:)
    integer :: i
    i = find(z, name)
    if (z%m(i)%ndim > 1) call fail(name//' is not one-dimensional')
    if (z%m(i)%is_int) then
      a = real(z%m(i)%k, dp)
    else
      a = z%m(i)%r
    end if
  end subroutine npz_real1

  subroutine npz_real2(z, name, a)
    type(npz_t), intent(in) :: z
    character(len=*), intent(in) :: name
    real(dp), allocatable, intent(out) :: a(:, :)
    integer :: i
    i = find(z, name)
    if (z%m(i)%ndim /= 2 .or. z%m(i)%is_int) call fail(name//' is not a 2-D float array')
    a = reshape(z%m(i)%r, [int(z%m(i)%shape(2)), int(z%m(i)%shape(1))])
  end subroutine npz_real2

  subroutine npz_int1(z, name, a)
    type(npz_t), intent(in) :: z
    character(len=*), intent(in) :: name
    integer, allocatable, intent(out) :: a(:)
    integer :: i
    i = find(z, name)
    if (.not. z%m(i)%is_int .or. z%m(i)%ndim > 1) call fail(name//' is not a 1-D int array')
    a = int(z%m(i)%k)
  end subroutine npz_int1

  real(dp) function npz_scalar(z, name) result(v)
    type(npz_t), intent(in) :: z
    character(len=*), intent(in) :: name
    real(dp), allocatable :: a(:)
    call npz_real1(z, name, a)
    v = a(1)
  end function npz_scalar

  pure integer(i8) function u16(b, p)
    character(len=*), intent(in) :: b
    integer(i8), intent(in) :: p
    u16 = int(ichar(b(p:p)), i8) + 256_i8 * int(ichar(b(p + 1:p + 1)), i8)
  end function u16

  pure integer(i8) function u32(b, p)
    character(len=*), intent(in) :: b
    integer(i8), intent(in) :: p
    u32 = u16(b, p) + 65536_i8 * u16(b, p + 2)
  end function u32

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    write(*, '(a)') 'npz_io: '//msg
    error stop 1
  end subroutine fail

end module npz_io
