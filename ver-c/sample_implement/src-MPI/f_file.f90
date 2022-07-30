! ************************************************************************** 
     subroutine file_open
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

   open (45,file='../data_out/force_pressure.dat' &
  &                                       ,status='unknown',form='formatted')
    write(45,*) &
   &  '       time, ID, Cx, Cy, Ct (Pressure), x_position, y_postion, angle'

   open (46,file='../data_out/force_friction.dat' &
  &                                       ,status='unknown',form='formatted')
    write(46,*) &
   &  '       time, ID, Cx, Cy, Ct (Friction), x_position, y_postion, angle'

   open (47,file='../data_out/force_total.dat' &
  &                                       ,status='unknown',form='formatted')
    write(47,*) &
   &  '       time, ID, Cx, Cy,  Ct (Total)  , x_position, y_postion, angle'

   open (48,file='../data_out/force_total_all.dat' &
  &                                       ,status='unknown',form='formatted')
    write(48,*) '       time, Cx, Cy,  Ct (Total)'

   open (55,file='../data_out/number_of_elements.dat' &
  &                                       ,status='unknown',form='formatted')
    write(55,*) '   time,    nvor_b,    nvor_s,   nvor_b + nvor_s '

   open (56,file='../data_out/check_elapsed_time.dat' &
  &                                       ,status='unknown',form='formatted')
    write(56,*) &
   & 'total_time,  matrix_v,  matrix_p,  get_velo,  core,  drift,  nascent, &
   &  trans '

  return
  end


! ************************************************************************** 
     subroutine file_write
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  implicit none

   write(55,'(f8.5,1x,i10,1x,i10,9x,i10)') &
  &     time, nvor_b, nvor_s, nvor_b+nvor_s

   write(56,'(e10.3, 2x, e10.3, 1x, e10.3, e10.3, e9.2, e9.2, e9.2, e9.2)') &
  &     total_time, time_matrix_velo, time_matrix_pres, time_get_velo,      &
  &     time_core,  time_drift,       time_nascent, time_trans

  return
  end


! ************************************************************************** 
     subroutine file_close
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

   close(45)
   close(46)
   close(47)
   close(48)
   close(55)
   close(56)

  return
  end


