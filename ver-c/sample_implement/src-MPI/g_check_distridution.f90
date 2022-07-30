! ************************************************************************** 
      subroutine check_distridution
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! ‰Q—v‘fˆÊ’u‚ÌŠm”F                                                           
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

!////////// blob elemnts //////////*

    new_b = 0

   do i = 1, nvor_b

    x0 = vor_b(1,i)
    y0 = vor_b(2,i)
    r0 = dsqrt(x0*x0 + y0*y0)

    cir2 = cir_b(i)
    cor2 = cor_b(i)

    call cross2(x1,y1,cir2,cor2,ic,height,icross)

    if (icross .eq. 1) cycle
    if (icross .eq. 0) then
     new_b = new_b + 1
     vor_b   (1,new_b) = x1
     vor_b   (2,new_b) = y1
     cir_b     (new_b) = cir2
     cor_b     (new_b) = cor2
     vor_vb(1,1,new_b) = vor_vb(1,1,i)
     vor_vb(2,1,new_b) = vor_vb(2,1,i)
     vor_vb(1,2,new_b) = vor_vb(1,2,i)
     vor_vb(2,2,new_b) = vor_vb(2,2,i)
     vor_vb(1,3,new_b) = vor_vb(1,3,i)
     vor_vb(2,3,new_b) = vor_vb(2,3,i)
     blob_id   (new_b) = blob_id   (i)
    end if

   end do

      nvor_b = new_b


!////////// sheet elemnts //////////*

    new_b = 0
    new_s = 0

    r_unit = 4.0d0*dh

   do i = 1, nvor_s

    x0 = vor_sc(1,i)
    y0 = vor_sc(2,i)
    r0 = dsqrt(x0*x0 + y0*y0)

    cor2 = cor_s(i)
    cir2 = cir_s(i)

!    vmx = 1.5d0*vor_vs(1,1,i) - 0.5d0*vor_vs(1,2,i)
!    vmy = 1.5d0*vor_vs(2,1,i) - 0.5d0*vor_vs(2,2,i)

    vmx=  23.0d0/12.0d0*vor_vs(1,1,i) &
   &     - 4.0d0/ 3.0d0*vor_vs(1,2,i) &
   &     + 5.0d0/12.0d0*vor_vs(1,3,i) 

    vmy=  23.0d0/12.0d0*vor_vs(2,1,i) &
   &     - 4.0d0/ 3.0d0*vor_vs(2,2,i) 
   &     + 5.0d0/12.0d0*vor_vs(2,3,i) 

    x1 = x0 + vmx*dt
    y1 = y0 + vmy*dt

    call cross2(x1,y1,cir2,cor2,ic,height,icross)

    if (icross .eq. 1) cycle

    if (icross .eq. 0) then

     c1 = 2.0d0*cor_s(i)
     c2 = vor_sr(i)

     if (height .le. r_unit .and. c1 .lt. c2) then

      new_s = new_s + 1

      gx = poi(1,1,ic) - poi(1,2,ic)
      gy = poi(2,1,ic) - poi(2,2,ic)
      gr=dsqrt(gx*gx + gy*gy)

      vor_sc  (1,new_s) = x1
      vor_sc  (2,new_s) = y1
      vor_s (1,1,new_s) = vor_sc(1,new_s) + 0.5d0*vor_sr(i)*gx/gr
      vor_s (2,1,new_s) = vor_sc(2,new_s) + 0.5d0*vor_sr(i)*gy/gr
      vor_s (1,2,new_s) = vor_sc(1,new_s) - 0.5d0*vor_sr(i)*gx/gr
      vor_s (2,2,new_s) = vor_sc(2,new_s) - 0.5d0*vor_sr(i)*gy/gr

      vor_vs(1,1,new_s) = vor_vs(1,1,i)
      vor_vs(2,1,new_s) = vor_vs(2,1,i)
      vor_vs(1,2,new_s) = vor_vs(1,2,i)
      vor_vs(2,2,new_s) = vor_vs(2,2,i)
      vor_vs(1,3,new_s) = vor_vs(1,3,i)
      vor_vs(2,3,new_s) = vor_vs(2,3,i)

      cir_s     (new_s) = cir2
      cor_s     (new_s) = cor2
      vor_sr    (new_s) = vor_sr    (i)
      sheet_id  (new_s) = sheet_id  (i)

     else

      new_b = new_b + 1

      vor_nb   (1,new_b) = x1
      vor_nb   (2,new_b) = y1
      vor_nvb(1,1,new_b) = vor_vs(1,1,i)
      vor_nvb(2,1,new_b) = vor_vs(2,1,i)
      vor_nvb(1,2,new_b) = vor_vs(1,2,i)
      vor_nvb(2,2,new_b) = vor_vs(2,2,i)
      vor_nvb(1,3,new_b) = vor_vs(1,3,i)
      vor_nvb(2,3,new_b) = vor_vs(2,3,i)
      cir_nb     (new_b) = cir_s     (i)
      cor_nb     (new_b) = dsqrt(2.0d0*vor_sr(i)*cor_s(i)/pi)
      blob_nid   (new_b) = sheet_id  (i)

     end if

    end if


   end do

      nvor_s = new_s


!////////// Recheck number of vortex elements //////////*

   do i = 1, new_b
    nvor_b = nvor_b + 1
    vor_b   (1,nvor_b) = vor_nb   (1,i)
    vor_b   (2,nvor_b) = vor_nb   (2,i)
    vor_vb(1,1,nvor_b) = vor_nvb(1,1,i)
    vor_vb(2,1,nvor_b) = vor_nvb(2,1,i)
    vor_vb(1,2,nvor_b) = vor_nvb(1,2,i)
    vor_vb(2,2,nvor_b) = vor_nvb(2,2,i)
    vor_vb(1,3,nvor_b) = vor_nvb(1,3,i)
    vor_vb(2,3,nvor_b) = vor_nvb(2,3,i)
    cir_b     (nvor_b) = cir_nb     (i)
    cor_b     (nvor_b) = cor_nb     (i)
    blob_id   (nvor_b) = blob_nid   (i)
   end do

!///////////////////////////////////*


   return
   end


