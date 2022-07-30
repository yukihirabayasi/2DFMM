
!************************
!-----// MODULES  //-----
!************************

!====================================================================
 module parameter                                                    
  implicit none                                                      
  integer, parameter :: n_id    = 7                                  
  integer, parameter :: n_lay   = 21                                 
  integer, parameter :: n_panel = 5001                               
  integer, parameter :: n_vor_b = 3000001                            
  integer, parameter :: n_vor_s = 10001                              
  integer, parameter :: n_sp    = 3000001                            
 end module parameter                                                
!====================================================================

!====================================================================
 module math_constant                                                
  implicit none                                                      
  double precision, parameter :: pi     = 3.141592653589793238       
  double precision, parameter :: pitwo  = 0.50d0 / pi                
  double precision, parameter :: pifour = 0.25d0 / pi                
  double precision, parameter :: g      = 9.80665d0                  
 end module math_constant                                            
!====================================================================

!====================================================================
 module gauss_legendre                                               
  ! --------------------------------------------------------- !      
  !  Abscissas and weights for the Gauss-Legendre quadrature  !      
  !  with "n_weight" points                                   !      
  ! --------------------------------------------------------- !      
   implicit none                                                     
                                                                     
! --- <<< This table contains values for n_weight = 8 >>> ---        
!   integer, parameter :: n_weight = 8                               
!   double precision, save :: phai (n_weight)                        
!   double precision, save :: w    (n_weight)                        
!  contains                                                          
!   subroutine intp_data                                             
!   phai(1) =  0.9602898564                                          
!   phai(2) =  0.7966664774                                          
!   phai(3) =  0.5255324099                                          
!   phai(4) =  0.1834346424                                          
!   phai(5) = -0.9602898564                                          
!   phai(6) = -0.7966664774                                          
!   phai(7) = -0.5255324099                                          
!   phai(8) = -0.1834346424                                          
!   w   (1) =  0.1012285362                                          
!   w   (2) =  0.2223810344                                          
!   w   (3) =  0.3137066458                                          
!   w   (4) =  0.3626837833                                          
!   w   (5) =  0.1012285362                                          
!   w   (6) =  0.2223810344                                          
!   w   (7) =  0.3137066458                                          
!   w   (8) =  0.3626837833                                          
!   end subroutine intp_data                                         
                                                                     
! --- <<< This table contains values for n_weight = 20 >>> ---       
   integer, parameter :: n_weight = 20                               
   double precision, save :: phai (n_weight)                         
   double precision, save :: w    (n_weight)                         
  contains                                                           
   subroutine intp_data                                              
   phai (1) =  0.993128599185094924786                               
   phai (2) =  0.963971927277913791268                               
   phai (3) =  0.912234428251325905868                               
   phai (4) =  0.839116971822218823395                               
   phai (5) =  0.746331906460150792614                               
   phai (6) =  0.636053680726515025453                               
   phai (7) =  0.510867001950827098004                               
   phai (8) =  0.373706088715419560673                               
   phai (9) =  0.227785851141645078080                               
   phai(10) =  0.076526521133497333755                               
   phai(11) = -0.993128599185094924786                               
   phai(12) = -0.963971927277913791268                               
   phai(13) = -0.912234428251325905868                               
   phai(14) = -0.839116971822218823395                               
   phai(15) = -0.746331906460150792614                               
   phai(16) = -0.636053680726515025453                               
   phai(17) = -0.510867001950827098004                               
   phai(18) = -0.373706088715419560673                               
   phai(19) = -0.227785851141645078080                               
   phai(20) = -0.076526521133497333755                               
   w    (1) =  0.017614007139152118312                               
   w    (2) =  0.040601429800386941331                               
   w    (3) =  0.062672048334109063570                               
   w    (4) =  0.083276741576704748725                               
   w    (5) =  0.101930119817240435037                               
   w    (6) =  0.118194531961518417312                               
   w    (7) =  0.131688638449176626898                               
   w    (8) =  0.142096109318382051329                               
   w    (9) =  0.149172986472603746788                               
   w   (10) =  0.152753387130725850698                               
   w   (11) =  0.017614007139152118312                               
   w   (12) =  0.040601429800386941331                               
   w   (13) =  0.062672048334109063570                               
   w   (14) =  0.083276741576704748725                               
   w   (15) =  0.101930119817240435037                               
   w   (16) =  0.118194531961518417312                               
   w   (17) =  0.131688638449176626898                               
   w   (18) =  0.142096109318382051329                               
   w   (19) =  0.149172986472603746788                               
   w   (20) =  0.152753387130725850698                               
   end subroutine intp_data                                          
                                                                     
 end module gauss_legendre                                           
!====================================================================

!====================================================================
 module global                                                       
                                                                     
  use parameter                                                      
  implicit none                                                      
                                                                     
  character(50), save :: file1                                       
  character(50), save :: file2                                       
  character(50), save :: file3                                       
  character(50), save :: file4                                       
                                                                     
! *** calculation parameters ***                                     
                                                                     
  integer, save :: nt                                                
  integer, save :: loop                                              
  integer, save :: i_write                                           
  integer, save :: i_continue                                        
  integer, save :: i_panel_height_check                              
  integer, save :: i_change_cond                                     
  integer, save :: i_kutta_cond                                      
  integer, save :: i_mirror                                          
  integer, save :: i_method                                          
  integer, save :: i_edge                                            
  integer, save :: i_layer                                           
                                                                     
  double precision, save :: d_mirror                                 
  double precision, save :: omz                                      
  double precision, save :: dt                                       
  double precision, save :: dt_base                                  
  double precision, save :: time                                     
  double precision, save :: re                                       
  double precision, save :: c_h                                      
  double precision, save :: ch_c                                     
  double precision, save :: c_dif                                    
  double precision, save :: cut_r                                    
  double precision, save :: uinf                                     
  double precision, save :: vinf                                     
  double precision, save :: theta                                    
  double precision, save :: theta0                                   
  double precision, save :: phi                                      
                                                                     
! *** parameter for wind turbine ***                                 
                                                                     
  double precision, save :: r_rot                                    
  double precision, save :: v_rot                                    
  double precision, save :: tsr                                      
                                                                     
! *** parameter for flapping ***                                     
                                                                     
  double precision, save :: cycle_time                               
  double precision, save :: pitch_angle_base_up                      
  double precision, save :: pitch_angle_base_down                    
  double precision, save :: dt_t                                     
  double precision, save :: dt_r                                     
  double precision, save :: dt_lag                                   
  double precision, save :: t_t                                      
  double precision, save :: t_r                                      
  double precision, save :: d_alpha_0                                
                                                                     
! *** parameter for move_aerofoil ***                                
                                                                     
  double precision, save :: velocity                                 
  double precision, save :: accelerate                               
  double precision, save :: velocity_base                            
  double precision, save :: uniform_velocity                         
  double precision, save :: u_accelerate                             
  double precision, save :: uniform_velocity_base                    
                                                                     
! *** parameter for move heaving ***                                 
                                                                     
  double precision, save :: amp                                      
  double precision, save :: freq                                     
  double precision, save :: phai_h                                   
                                                                     
! *** parameter for move feathering ***                              
                                                                     
  double precision, save :: s                                        
  double precision, save :: s_int                                    
  double precision, save :: s_amp                                    
  double precision, save :: s_old                                    
  double precision, save :: k_p                                      
  double precision, save :: omega_f                                  
  double precision, save :: ds_p                                     
  double precision, save :: f_cx                                     
  double precision, save :: f_cy                                     
								     
! *** parts variables ***                                            
                                                                     
  integer, save :: edge (2,n_id)                                     
                                                                     
  double precision, save :: edge_cir       (2,n_id)                  
  double precision, save :: cir_r            (n_id)                  
  double precision, save :: phase            (n_id)                  
  double precision, save :: pgw           (3, n_id)                  
  double precision, save :: parts_position(3, n_id)                  
                                                                     
! *** panel element ***                                              
                                                                     
  integer, save :: nw                                                
  integer, save :: npanel                                            
  integer, save :: panel_id   (n_panel)                              
  integer, save :: panel_n  (2,n_panel)                              
                                                                     
  double precision, save :: q             (n_panel)                  
  double precision, save :: ph            (n_panel)                  
  double precision, save :: poi     (2, 2, n_panel)                  
  double precision, save :: poi_b   (2, 2, n_panel)                  
  double precision, save :: uw      (2,    n_panel)                  
  double precision, save :: cir_g         (n_panel)                  
  double precision, save :: poi_c      (2, n_panel)                  
  double precision, save :: ds            (n_panel)                  
  double precision, save :: poi_n      (2, n_panel)                  
  double precision, save :: cp_x          (n_panel)                  
  double precision, save :: cp_y          (n_panel)                  
  double precision, save :: cf_x          (n_panel)                  
  double precision, save :: cf_y          (n_panel)                  
  double precision, save :: cp_t          (n_panel)                  
  double precision, save :: cf            (n_panel)                  
  double precision, save :: vn_drft       (n_panel)                  
  double precision, save :: gam_drft      (n_panel)                  
  double precision, save :: rot1          (n_panel)                  
  double precision, save :: rot2          (n_panel)                  
  double precision, save :: rot3          (n_panel)                  
  double precision, save :: rot4          (n_panel)                  
  double precision, save :: xc      (2, 2, n_panel)                  
  double precision, save :: vm      (2, 2, n_panel)                  
  double precision, save :: am      (2, 2, n_panel)                  
  double precision, save :: dh                                       
  double precision, save :: vn_vis                                   
  double precision, save :: ds_t                                     
                                                                     
! *** Vortex element (blob) ***                                      
                                                                     
  integer, save :: nvor_b                                            
  integer, save :: blob_id  (n_vor_b)                                
  integer, save :: blob_nid (n_vor_b)                                
                                                                     
  double precision, save :: cir_b         (n_vor_b)                  
  double precision, save :: cor_b         (n_vor_b)                  
  double precision, save :: vor_b      (2, n_vor_b)                  
  double precision, save :: vor_vb  (2, 3, n_vor_b)                  
  double precision, save :: cir_nb        (n_vor_b)                  
  double precision, save :: cor_nb        (n_vor_b)                  
  double precision, save :: vor_nb     (2, n_vor_b)                  
  double precision, save :: vor_nvb (2, 3, n_vor_b)                  
                                                                     
! *** layer ***                                                      
                                                                     
  integer, save :: nlayer                    (n_panel)               
  integer, save :: lay_id             (n_lay, n_panel)               
                                                                     
  double precision, save :: ypls             (n_panel)               
  double precision, save :: lay_h            (n_panel)               
  double precision, save :: lay_x  (2, n_lay, n_panel)               
  double precision, save :: lay_y  (2, n_lay, n_panel)               
  double precision, save :: lay_vx (2, n_lay, n_panel)               
  double precision, save :: lay_vy (2, n_lay, n_panel)               
  double precision, save :: lay_c  (2, n_lay, n_panel)               
  double precision, save :: lay_vc (2, n_lay, n_panel)               
                                                                     
! *** Vortex element (sheet) ***                                     
                                                                     
  integer, save :: nvor_s                                            
  integer, save :: sheet_id      (n_vor_s)                           
                                                                     
  double precision, save :: vor_s   (2, 2, n_vor_s)                  
  double precision, save :: vor_sc     (2, n_vor_s)                  
  double precision, save :: vor_vs  (2, 3, n_vor_s)                  
  double precision, save :: vor_sr        (n_vor_s)                  
  double precision, save :: cir_s         (n_vor_s)                  
  double precision, save :: cor_s         (n_vor_s)                  
                                                                     
! *** elapsed time check ***                                         
                                                                     
  integer,save :: hz                                                 
  integer,save :: clock_1                                            
  integer,save :: clock_2                                            
                                                                     
  double precision,save :: total_time                                
  double precision,save :: time_matrix_velo                          
  double precision,save :: time_matrix_pres                          
  double precision,save :: time_get_velo                             
  double precision,save :: time_core                                 
  double precision,save :: time_drift                                
  double precision,save :: time_nascent                              
  double precision,save :: time_trans                                
                                                                     
! *** parallel computation ***                                       
                                                                     
  integer, save :: my_rank                                           
  integer, save :: num_procs                                         
  integer, save :: ierr                                              
  integer, save :: fquit                                             
                                                                     
 end module global                                                   
!====================================================================


