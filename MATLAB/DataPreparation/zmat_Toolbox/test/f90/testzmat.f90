!------------------------------------------------------------------------------
!  ZMatLib test program
!------------------------------------------------------------------------------
!
!> @author Qianqian Fang <q.fang at neu.edu>
!> @brief A demo program to call zmatlib functions to encode and decode
!
! DESCRIPTION:
!> This demo program shows how to call zmat_run/zmat_encode/zmat_decode to
!> perform buffer encoding and decoding
!------------------------------------------------------------------------------

program zmatdemo
!------------------------------------------------------------------------------
! step 1: add the below line to use zmatlib unit
!------------------------------------------------------------------------------
use zmatlib
use iso_c_binding, only: c_double
implicit none

!------------------------------------------------------------------------------
! step 2: define needed input and output variables
!------------------------------------------------------------------------------

! here are the fortran arrays/variables that you want to pass to zmat_run
character (len=128), target :: inputstr
real(kind=c_double), target :: dat(10)

! here are the interface variables, note that inputbuf/outputbuf are pointers (void*)
integer(kind=c_int) :: ret, res, i
integer(kind=c_size_t) :: inputlen, outputlen
type(c_ptr) :: inputbuf, outputbuf

! zmat output is defined as a char pointer, one must free this after each call
character(kind=c_char),pointer :: myout(:)

!------------------------------------------------------------------------------
! step 3: set whatever fortran variables
!------------------------------------------------------------------------------
inputstr="__o000o__(o)(o)__o000o__ =^_^=  __o000o__(o)(o)__o000o__"

print *, trim(inputstr)

!------------------------------------------------------------------------------
! step 4: set inputbuf to the variable you want to encode, then call zmat_run
!------------------------------------------------------------------------------

inputbuf=c_loc(inputstr)
inputlen=len(trim(inputstr))

res=zmat_run(inputlen,inputbuf,outputlen, outputbuf, zmBase64, ret, 1);

!------------------------------------------------------------------------------
! step 5: once done, call c_f_pointer to transfer myout to outputbuf, print results
!------------------------------------------------------------------------------

call c_f_pointer(outputbuf, myout, [outputlen])
print *, myout

!------------------------------------------------------------------------------
! step 6: free the c-allocated buffer
!------------------------------------------------------------------------------

call zmat_free(outputbuf)

!------------------------------------------------------------------------------
! repeat step 3 for another variable: set whatever fortran variables
!------------------------------------------------------------------------------
dat = (/ (i, i = 1, 10) /)

!------------------------------------------------------------------------------
! repeat step 4: set inputbuf to the variable you want to encode, then call zmat_run
!------------------------------------------------------------------------------

inputlen=sizeof(dat)
inputbuf=c_loc(dat)

res=zmat_run(inputlen,inputbuf,outputlen, outputbuf, zmBase64, ret, 1);

!------------------------------------------------------------------------------
! repeat step 5: once done, call c_f_pointer to transfer myout to outputbuf, print results
!------------------------------------------------------------------------------

call c_f_pointer(outputbuf, myout, [outputlen])
print *, myout

!------------------------------------------------------------------------------
! repeat step 6: free the c-allocated buffer
!------------------------------------------------------------------------------

call zmat_free(outputbuf)

end program
