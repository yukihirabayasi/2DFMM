! ************************************************************************** 
      subroutine get_velo
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! Calculation of convective velocity of vortex elements.                     
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

   do i = 1, nvor_b
    x1 = vor_b(1,i)
    y1 = vor_b(2,i)
    call indus_s(x1,y1,u1,v1)
    call indus_v(x1,y1,u2,v2)
    call biot_b (x1,y1,u3,v3)
    call biot_s (x1,y1,u4,v4)
    vmx = uinf + (-omz*y1)
    vmy = vinf + ( omz*x1)
    vor_vb(1,1,i) = vmx + u1 + u2 + u3 + u4
    vor_vb(2,1,i) = vmy + v1 + v2 + v3 + v4
   end do

   do i = 1, nvor_s
    x1 = vor_sc(1,i)
    y1 = vor_sc(2,i)
    call indus_s(x1,y1,u1,v1)
    call indus_v(x1,y1,u2,v2)
    call biot_b (x1,y1,u3,v3)
    call biot_s (x1,y1,u4,v4)
    vmx = uinf + (-omz*y1)
    vmy = vinf + ( omz*x1)
    vor_vs(1,1,i) = vmx + u1 + u2 + u3 + u4
    vor_vs(2,1,i) = vmy + v1 + v2 + v3 + v4
   end do


   return
   end


