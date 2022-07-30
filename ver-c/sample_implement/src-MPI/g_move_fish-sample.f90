! ************************************************************************** 
      subroutine move
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

   double precision :: a0(5)
   double precision :: tm(5)
   double precision :: fp(5)
   double precision :: a1

   data   a0/ 6.570d0, 0.026d0, 0.026d0, 0.026d0, 0.026d0 /
   data   tm/ 0.420d0, 0.650d0, 0.650d0, 0.650d0, 0.650d0 /
   data   fp/ 50.00d0, 0.650d0, 0.650d0, 0.650d0, 0.650d0 /

!/////////// äÓñ{óVâjèåè ///////////*
!----- a0 (êUïùäp[6.570])       -----*
!----- tm (ñ≥éüå≥é¸ä˙[0.420])   -----*
!----- fp (îˆïîà ëäíxÇÍ[50.00]) -----*
!////////////////////////////////////*


!------------------
  do i = 1, npanel 
!------------------

   k = panel_id(i)

   fre = 1.0d0/tm(k)
   fai_p = fp(k)*pi/180.0d0

!-------------
  do j = 1, 2 
!-------------

   xc(1,j,i) = poi_b(1,j,i) - pgw(1,k)
   x  = xc(1,j,i)
   r  = dabs(x)
   a1 = a0(k)*pi/180.0d0

   if (x .le. 0.0) then
    theta_p = a1*dcos(2.0d0*pi*fre*time - phase(k) + pi)
   else
    theta_p = a1*dcos(2.0d0*pi*fre*time - phase(k)-(r/0.75d0)*fai_p)
   end if

    dx = r*dcos(theta_p)
    dy = r*dsin(theta_p)

    poi  (1,j,i) = dx
    poi  (2,j,i) = poi_b(2,j,i) +dy

!--------------
   end do
!--------------

!////////// ë¨ìxÇ∆â¡ë¨ìxÇÃåvéZ //////////*

    x1 = poi_b(1,1,i) - pgw(1,k)
    x2 = poi_b(1,2,i) - pgw(1,k)
    x  = 0.5*(x1+x2)
    r  = dabs(x)

    vm(1,2,i) = vm(1,1,i)
    vm(2,2,i) = vm(2,1,i)

    am(1,2,i) = am(1,1,i)
    am(2,2,i) = am(2,1,i)

    if (x .le. 0.0d0) then
     theta_p = a1*dcos(2.0d0*pi*fre*time - phase(k)) + pi
     omega_p = -2.0d0*pi*fre*a1*dsin(2.0d0*pi*fre*time - phase(k))
     alpha_p = -(2.0d0*pi*fre)**2*a1*dcos(2.0d0*pi*fre*time - phase(k))
    else
     theta_p = a1*dcos(2.0d0*pi*fre*time - phase(k) &
               -(r/0.75d0)*fai_p)
     omega_p = -2.0d0*pi*fre*a1*dsin(2.0d0*pi*fre*time - phase(k) &
               -(r/0.75d0)*fai_p)
     alpha_p = -(2.0d0*pi*fre)**2*a1*dcos(2.0d0*pi*fre*time - phase(k) &
               -(r/0.75d0)*fai_p)
    end if

     vm(1,1,i) = -omega_p*r*dsin(theta_p)
     vm(2,1,i) =  omega_p*r*dcos(theta_p)

     am(1,1,i) = -alpha_p*r*dsin(theta_p) - (omega_p*omega_p)*r*dcos(theta_p)
     am(2,1,i) =  alpha_p*r*dcos(theta_p) - (omega_p*omega_p)*r*dsin(theta_p)

!---------
   end do 
!---------

   return
   end

