module MODE_WRITE

  contains
  subroutine WRITE_2D(var,fileid)
    ! 
    ! BC created August 2021
    !  Write 2D variable var into fileid (append mode)
    ! @TODO: interface to handle non-real variables.
    implicit none
    real, dimension(:,:), intent(in):: var
    integer, intent(in) :: fileid
    integer :: pos,i,j
    inquire(unit = fileid, pos = pos)
    write(fileid, pos = pos) ((var(i,j),j=1,size(var,2)),i=1,size(var,1))
  end subroutine WRITE_2D

end module MODE_WRITE

module IOUNITS
  integer, parameter :: &
    udmp = 11,         &! Dump file unit number
    umet = 21,         &! Driving file unit number
    uout = 31,         &! Output file unit number
    ustr = 41           ! Start file unit number
end module IOUNITS

module FILES
  character(len=200) :: &
    dump_file           ! Dump file name
end module FILES
