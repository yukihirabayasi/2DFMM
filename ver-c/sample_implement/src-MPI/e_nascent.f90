! ************************************************************************** 
     subroutine nascent
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

   double precision :: coef_omg
   double precision :: gam_vis

!  ////////// calculate flux ///////////

   c_dif  = 1.136d0
   vn_vis = 0.5d0 * c_dif**2 / dh / re

   call flux_drft

!  ////////// nascent vortex elements ///////////

   do i = 1, npanel

     v_flux = vn_drft(i) + vn_vis

     if ( v_flux .lt. 0.0d0 ) cycle

!  ////////// circulation of viscosity effect ///////////

     ! j = 1
       ux = lay_vc(1,1,i)
       uy = lay_vc(2,1,i)
       vtt = ( uy*poi_n(1,i) - ux*poi_n(2,i) )

     if ( i_layer .eq. 1 ) then
       do j = 2, nlayer(i)
         ux = lay_vc(1,j,i) - lay_vc(1,j-1,i)
         uy = lay_vc(2,j,i) - lay_vc(2,j-1,i)
         vtt = vtt + ( uy*poi_n(1,i) - ux*poi_n(2,i) )
       end do
     end if

     omg_z = vtt / dh

!  ////////// profile of nascent vortex element ///////////

     dh2 = v_flux * dt

     coef_omg = dh / ( dh + vn_vis * dt )
     gam_vis  = ( coef_omg * omg_z )*( vn_vis * ds(i) * dt )

     nvor_s = nvor_s + 1

     cir_s (nvor_s) = gam_vis + gam_drft(i)
     cor_s (nvor_s) = 0.5d0 * dh2
     vor_sr(nvor_s) = ds(i)

     do j = 1, nw
       if ( edge(1,j) .eq. i ) edge_cir(1,j) = cir_s(nvor_s)
       if ( edge(2,j) .eq. i ) edge_cir(2,j) = cir_s(nvor_s)
     end do

     dh1 = dh + 2.0d0*cor_s(nvor_s)

     vor_s (1,1,nvor_s) =  poi(1,1,i) + dh1 * poi_n(1,i)
     vor_s (2,1,nvor_s) =  poi(2,1,i) + dh1 * poi_n(2,i)
     vor_s (1,2,nvor_s) =  poi(1,2,i) + dh1 * poi_n(1,i)
     vor_s (2,2,nvor_s) =  poi(2,2,i) + dh1 * poi_n(2,i)

     vor_sc  (1,nvor_s) = 0.5d0 * ( vor_s(1,1,nvor_s) + vor_s(1,2,nvor_s) )
     vor_sc  (2,nvor_s) = 0.5d0 * ( vor_s(2,1,nvor_s) + vor_s(2,2,nvor_s) )

     vor_vs(1,1,nvor_s) = lay_vc(1,nlayer(i),i) + vm(1,1,i)
     vor_vs(2,1,nvor_s) = lay_vc(2,nlayer(i),i) + vm(2,1,i)

     vor_vs(1,2,nvor_s) = vor_vs(1,1,nvor_s)
     vor_vs(2,2,nvor_s) = vor_vs(2,1,nvor_s)

     vor_vs(1,3,nvor_s) = vor_vs(1,1,nvor_s)
     vor_vs(2,3,nvor_s) = vor_vs(2,1,nvor_s)

     sheet_id  (nvor_s) = panel_id(i)

   end do

  return
  end


! ************************************************************************** 
     subroutine flux_drft
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

   double precision :: gam_c (2)
   double precision :: sgam_t(2)

   do i = 1, npanel

     do k = 1, 2

       ! j = 1
         rx = lay_x(k,1,i) - lay_c(1,1,i)
         ry = lay_y(k,1,i) - lay_c(2,1,i)
         r  = dsqrt( rx**2 + ry**2 )

         ux  = lay_vx(k,1,i)
         uy  = lay_vy(k,1,i)

         sgam_t(k) = ( ux*rx + uy*ry ) / r * 0.5d0 * lay_h(i)
         gam_c (k) = (( ux*rx + uy*ry ) / r )**2 * 0.5d0	! Caution

       if ( i_layer .eq. 1 ) then
         do j = 2, nlayer(i)
           rx = lay_x(k,j,i) - lay_c(1,j,i)
           ry = lay_y(k,j,i) - lay_c(2,j,i)
           r  = dsqrt( rx**2 + ry**2 )

           ux1  = lay_vx(k,j,i)
           uy1  = lay_vy(k,j,i)
           ux2  = lay_vx(k,j-1,i)
           uy2  = lay_vy(k,j-1,i)

           sgam_t(k) = sgam_t(k) + ( (ux1+ux2)*rx + (uy1+uy2)*ry ) / r * 0.5d0 * lay_h(i)
           gam_c (k) = gam_c (k) + (( (ux1-ux2)*rx + (uy1-uy2)*ry ) / r )**2 * 0.5d0	! Caution
         end do
       end if

     end do

     vn_drft (i) = -( sgam_t(1) + sgam_t(2) ) / ds(i)
     gam_drft(i) =  ( gam_c(1) - gam_c(2) ) * dt

   end do

  return
  end


! ************************************************************************** 
      subroutine kelvine
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! CHECK THE STRENGTH OF VORTEX CIRCULATION ( CIRCULATION THEOREM )           
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   sgam = 0.0d0
   ds_t = 0.0d0

   do i = 1, nvor_s
     sgam = sgam + cir_s(i)
   end do

   do i = 1, nvor_b
     sgam = sgam + cir_b(i)
   end do

   do i = 1, nw
     sgam = sgam + cir_r(i)
   end do

   do i = 1, npanel
     ds_t = ds_t + ds(i)
   end do

   do i = 1, npanel
     cir_g(i) = sgam / ds_t
   end do

  return
  end


