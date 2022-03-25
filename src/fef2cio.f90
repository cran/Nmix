  module fef2cio

! front-end to Peter Soerensen's f2cio module

  contains

!---------------------------------------------------

  subroutine openabf2cio(ifn,lfn,fp,append)

  use iso_c_binding
  use f2cio
  
  implicit none
  integer(c_int) :: lfn,ifn(lfn)  
  character(len=lfn, kind=c_char) :: fn
  character(len=20, kind=c_char) :: mode
  character(len=1000, kind=c_char) :: filename
  type(c_ptr):: fp
  integer(c_int):: i,append
  
  do i=1,lfn
    fn(i:i) = char(ifn(i))
  enddo

  filename = fn(1:lfn) // C_NULL_CHAR
  if(append.eq.1) then
	mode = 'ab' // C_NULL_CHAR
  else
	mode = 'wb' // C_NULL_CHAR
  end if
  fp = fopen(filename, mode)

  end subroutine openabf2cio

!---------------------------------------------------

  subroutine wdf2cio(x,n,fp)

  use iso_c_binding
  use f2cio
  
  implicit none
  type(c_ptr):: fp
  integer(c_int):: cfres,n

  integer, parameter :: dbl = kind(1.0d0)
  real(kind=dbl), target :: x(n)
  
  cfres = fwrite(c_loc(x),8,n,fp)

  end subroutine wdf2cio

!---------------------------------------------------

  subroutine wrf2cio(x,n,fp)

  use iso_c_binding
  use f2cio
  
  implicit none
  type(c_ptr):: fp
  integer(c_int):: cfres,n

  integer, parameter :: sng = kind(1.0)
  real(kind=sng), target :: x(n)
  
  cfres = fwrite(c_loc(x),4,n,fp)

  end subroutine wrf2cio

!---------------------------------------------------

  subroutine wif2cio(x,n,fp)

  use iso_c_binding
  use f2cio
  
  implicit none
  type(c_ptr):: fp
  integer(c_int):: cfres,n

  integer, target :: x(n)
  
  cfres = fwrite(c_loc(x),4,n,fp)

  end subroutine wif2cio

!---------------------------------------------------

  subroutine closef2cio(fp)

  use iso_c_binding
  use f2cio
  
  implicit none
  type(c_ptr):: fp
  integer(c_int):: cfres

  cfres = fclose(fp)

  end subroutine closef2cio

!---------------------------------------------------
     
  end module fef2cio

