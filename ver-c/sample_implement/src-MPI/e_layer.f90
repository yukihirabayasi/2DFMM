! ************************************************************************** 
     subroutine layer
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   character(40) :: nfile
   character(40) :: fout_h
   character(40) :: fout_t
   double precision :: laynx, layny
   double precision :: omg_l
   double precision :: uw_new
   double precision :: omg_l_1
   double precision :: ypls_1(npanel)

   fout_h = '../data_out/yplus_'
   fout_t = '_step.dat'

   write(nfile,'(A18,i5.5,A9)') fout_h, nt, fout_t

   open(36,file=nfile,status='unknown',form='formatted')

    write(36,*) '  id,    y+,          layer,       ux_1,        yplus_1'

    do i = 1, npanel

     omg_l   = dabs ( uw(2,i)*poi_n(1,i) - uw(1,i)*poi_n(2,i) ) / dh
     ypls(i) = dsqrt( omg_l*re ) * dh

     nlayer(i) = dint( ypls(i) ) + 1
     lay_h (i) = dh / nlayer(i)

     if ( i_layer .eq. 0 ) then
       nlayer(i) = 1
       lay_h (i) = dh
     end if

     do j = 1, nlayer(i)

       laynx = lay_h(i) * poi_n(1,i)
       layny = lay_h(i) * poi_n(2,i)

       if ( j .eq. 1 ) then
         lay_x(1,j,i) = poi(1,1,i) + laynx
         lay_y(1,j,i) = poi(2,1,i) + layny
         lay_x(2,j,i) = poi(1,2,i) + laynx
         lay_y(2,j,i) = poi(2,2,i) + layny
       else
         lay_x(1,j,i) = lay_x(1,j-1,i) + laynx
         lay_y(1,j,i) = lay_y(1,j-1,i) + layny
         lay_x(2,j,i) = lay_x(2,j-1,i) + laynx
         lay_y(2,j,i) = lay_y(2,j-1,i) + layny
       end if

       lay_c(1,j,i) = 0.5d0 * ( lay_x(1,j,i) + lay_x(2,j,i) )
       lay_c(2,j,i) = 0.5d0 * ( lay_y(1,j,i) + lay_y(2,j,i) )

       do k = 1, 2

         x1 = lay_x(k,j,i)
         y1 = lay_y(k,j,i)

         call indus_s( x1,y1,u1,v1 )
         call indus_v( x1,y1,u2,v2 )
         call biot_b ( x1,y1,u3,v3 )
         call biot_s ( x1,y1,u4,v4 )

         vmx = uinf + ( -omz*y1 ) - vm(1,1,i)
         vmy = vinf + (  omz*x1 ) - vm(2,1,i)
         lay_vx(k,j,i) = vmx + ( u1 + u2 + u3 + u4 )
         lay_vy(k,j,i) = vmy + ( v1 + v2 + v3 + v4 )

       end do

       lay_vc(1,j,i) = 0.5d0 * ( lay_vx(1,j,i) + lay_vx(2,j,i) )
       lay_vc(2,j,i) = 0.5d0 * ( lay_vy(1,j,i) + lay_vy(2,j,i) )
       lay_id  (j,i) = panel_id (i)

     end do

     omg_l_1   = dabs( lay_vc(2,1,i)*poi_n(1,i) - lay_vc(1,1,i)*poi_n(2,i) ) / lay_h(i)
     ypls_1(i) = dsqrt( omg_l_1 * re ) * lay_h(i)
     uw_new    = dsqrt( lay_vc(1,nlayer(i),i)**2 + lay_vc(1,nlayer(i),i)**2 )

     write(36,'(i5,30f13.6)') i, ypls(i), dble(nlayer(i)), uw_new, ypls_1(i)

   end do

  return
  end


