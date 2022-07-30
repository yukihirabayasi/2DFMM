! ************************************************************************** 
     subroutine kutta_condition
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! NASCENT OF VORTEX ELEMNTS ON TRAILING EDGE                                 
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

!===============
   do i = 1, nw 
!===============

   n1 = edge(1,i)
   n2 = edge(2,i)

   gam_u = lay_vy(1,nlayer(i),n1)*poi_n(1,n1) - lay_vx(1,nlayer(i),n1)*poi_n(2,n1)
   gam_d = lay_vy(2,nlayer(i),n2)*poi_n(1,n2) - lay_vx(2,nlayer(i),n2)*poi_n(2,n2)

   rx = poi(1,1,n1) - 0.5d0*( poi_c(1,n1) + poi_c(1,n2) )
   ry = poi(2,1,n1) - 0.5d0*( poi_c(2,n1) + poi_c(2,n2) )
   r  = dsqrt( rx**2 + ry**2 )

   u1 = ( lay_vx(1,nlayer(i),n1)*rx + lay_vy(1,nlayer(i),n1)*ry ) / r
   u2 = ( lay_vx(2,nlayer(i),n2)*rx + lay_vy(2,nlayer(i),n2)*ry ) / r

   px0 = 0.5d0 * ( poi(1,1,n1) + poi(1,2,n2) )
   py0 = 0.5d0 * ( poi(2,1,n2) + poi(2,2,n2) )

   dl = 0.5d0*( u1 + u2 )*dt

   if ( dl .gt. 1.0d-4 ) then
     nvor_b = nvor_b + 1
     cir_b     (nvor_b) = ( gam_u + gam_d )*dl - ( edge_cir(1,i) + edge_cir(2,i) )
     cor_b     (nvor_b) = 0.5d0 * dl
     vor_b   (1,nvor_b) = poi(1,1,n1) + cor_b(nvor_b) * rx / r
     vor_b   (2,nvor_b) = poi(2,1,n1) + cor_b(nvor_b) * ry / r
     vor_vb(1,1,nvor_b) = 0.5d0 * ( u1 + u2 ) * rx / r
     vor_vb(2,1,nvor_b) = 0.5d0 * ( u1 + u2 ) * ry / r
     vor_vb(1,2,nvor_b) = vor_vb(1,1,nvor_b)
     vor_vb(2,2,nvor_b) = vor_vb(2,1,nvor_b)
     vor_vb(1,3,nvor_b) = vor_vb(1,1,nvor_b)
     vor_vb(2,3,nvor_b) = vor_vb(2,1,nvor_b)
     blob_id   (nvor_b) = i
   end if

!   if ( dl .gt. 1.0d-4 ) then
!     nvor_s = nvor_s + 1
!     cir_s     (nvor_s) = ( gam_u + gam_d )*dl - ( edge_cir(1,i) + edge_cir(2,i) )
!     cor_s     (nvor_s) = 0.5d0 * dl
!     vor_sr    (nvor_s) = dl
!     vor_s (1,1,nvor_s) = poi(1,1,n1) +   cor_s(nvor_s)        * rx / r
!     vor_s (2,1,nvor_s) = poi(2,1,n1) +   cor_s(nvor_s)        * ry / r
!     vor_s (1,2,nvor_s) = poi(1,1,n1) + ( cor_s(nvor_s) + dl ) * rx / r
!     vor_s (2,2,nvor_s) = poi(2,1,n1) + ( cor_s(nvor_s) + dl ) * ry / r
!     vor_sc  (1,nvor_s) = 0.5d0 * ( vor_s(1,1,nvor_s) + vor_s(1,2,nvor_s) )
!     vor_sc  (2,nvor_s) = 0.5d0 * ( vor_s(2,1,nvor_s) + vor_s(2,2,nvor_s) )
!     vor_vs(1,1,nvor_s) = 0.5d0 * ( u1 + u2 ) * rx / r
!     vor_vs(2,1,nvor_s) = 0.5d0 * ( u1 + u2 ) * ry / r
!     vor_vs(1,2,nvor_s) = vor_vs(1,1,nvor_s)
!     vor_vs(2,2,nvor_s) = vor_vs(2,1,nvor_s)
!     vor_vs(1,3,nvor_s) = vor_vs(1,1,nvor_s)
!     vor_vs(2,3,nvor_s) = vor_vs(2,1,nvor_s)
!     sheet_id  (nvor_s) = i
!   end if

!=========
   end do 
!=========

  return
  end


