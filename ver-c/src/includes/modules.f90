
!************************
!-----// MODULES  //-----
!************************

!====================================================================
 module parameter                                                    
  implicit none                                                      
  integer, parameter :: n_id    = 7                                  
  integer, parameter :: n_lay   = 11                                 
  integer, parameter :: n_field = 2001                               
  integer, parameter :: n_panel = 5001                               
  integer, parameter :: n_vor_b = 3000001                            
  integer, parameter :: n_vor_s = 10001                              
  integer, parameter :: n_cell  = n_vor_b                            
  integer, parameter :: n_node  = n_vor_b*2                          
 end module parameter                                                
!====================================================================

!====================================================================
 module math_constant                                                
  implicit none                                                      
  real(8), parameter :: pi     = 3.141592653589793238                
  real(8), parameter :: pitwo  = 0.50d0 / pi                         
  real(8), parameter :: pifour = 0.25d0 / pi             
  double precision,save :: gpx                    
 end module math_constant                                            
!====================================================================

!====================================================================
 module gauss_legendre                                               
  ! --------------------------------------------------------- !      
  !  Abscissas and weights for the Gauss-Legendre quadrature  !      
  !  with "n_weight" points                                   !      
  ! --------------------------------------------------------- !      
   implicit none                                                     
                                                                     
! --- <<< This table contains values for n_weight=20 >>> ---         
   integer, parameter :: n_weight = 20                               
   real(8), save :: phai (n_weight)                                  
   real(8), save :: w    (n_weight)                                  
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
                                                                     
  character(40), save :: file_panel                                  
  character(40), save :: file_continue                               
                                                                     
! *** calculation parameters ***                                     
                                                                     
  integer, save :: nt                                                
  integer, save :: loop                                              
  integer, save :: i_write                                           
  integer, save :: i_continue                                        
  integer, save :: i_panel_height_check                              
  integer, save :: i_change_cond                                     
  integer, save :: i_kutta_cond                                      
  integer, save :: i_mirror                                          
  integer, save :: i_edge                                            
                                                                     
  real(8), save :: ch                                                
  real(8), save :: cut_r                                             
  real(8), save :: d_mirror                                          
  real(8), save :: dt                                                
  real(8), save :: dt_base                                           
  real(8), save :: omz                                               
  real(8), save :: pheta                                             
  real(8), save :: theta                                             
  real(8), save :: re                                                
  real(8), save :: time                                              
  real(8), save :: uinf                                              
  real(8), save :: vinf                                              
                                                                     
! *** parameter for wind turbine ***                                 
                                                                     
  real(8), save :: tsr                                               
  real(8), save :: r_rot                                             
  real(8), save :: v_rot                                             
                                                                     
! *** parameter for flapping ***                                     
                                                                     
  real(8), save :: cycle_time                                        
  real(8), save :: d_alpha_0                                         
  real(8), save :: dt_t                                              
  real(8), save :: dt_r                                              
  real(8), save :: dt_lag                                            
  real(8), save :: pitch_angle_base_up                               
  real(8), save :: pitch_angle_base_down                             
  real(8), save :: t_t                                               
  real(8), save :: t_r                                               
                                                                     
! *** parameter for move_aerofoil ***                                
                                                                     
  real(8), save :: accelerate                                        
  real(8), save :: velocity                                          
  real(8), save :: velocity_base                                     
  real(8), save :: u_accelerate                                      
  real(8), save :: uniform_velocity                                  
  real(8), save :: uniform_velocity_base                             
                                                                     
! *** parameter for move heaving ***                                 
                                                                     
  real(8), save :: amp                                               
  real(8), save :: freq                                              
                                                                     
! *** parameter for move feathering ***                              
                                                                     
  real(8), save :: k_p                                               
  real(8), save :: s_int                                             
  real(8), save :: s_amp                                             
  real(8), save :: f_cx                                              
  real(8), save :: f_cy                                              
                                                                     
! *** parameter for move rotating ***                                
                                                                     
  real(8), save :: omega_t                                           
                                                                     
! *** parts variables ***                                            
                                                                     
  integer, save :: edge                 (2, n_id)                    
                                                                     
  real(8), save :: cir_r                   (n_id)                    
  real(8), save :: edge_cir             (2, n_id)                    
  real(8), save :: parts_position       (3, n_id)                    
                                                                     
! *** panel element ***                                              
                                                                     
  integer, save :: nw                                                
  integer, save :: npanel                                            
  integer, save :: panel_id             (n_panel)                    
  integer, save :: panel_n           (2, n_panel)                    
                                                                     
  real(8), save :: q                    (n_panel)                    
  real(8), save :: ph                   (n_panel)                    
  real(8), save :: ds                   (n_panel)                    
  real(8), save :: cir_g                (n_panel)                    
  real(8), save :: uw                (2, n_panel)                    
  real(8), save :: poi            (2, 2, n_panel)                    
  real(8), save :: poi_b          (2, 2, n_panel)                    
  real(8), save :: poi_c             (2, n_panel)                    
  real(8), save :: poi_n             (2, n_panel)                    
  real(8), save :: cf                   (n_panel)                    
  real(8), save :: cp                   (n_panel)                    
  real(8), save :: vm             (2, 2, n_panel)                    
  real(8), save :: am             (2, 2, n_panel)                    
  real(8), save :: dh                                                
  real(8), save :: ds_t                                              
                                                                     
! *** multi layer ***                                                
                                                                     
  integer, save :: nlayer               (n_panel)                    
                                                                     
  real(8), save :: y_plus               (n_panel)                    
  real(8), save :: gam           (n_lay, n_panel)                    
  real(8), save :: vn            (n_lay, n_panel)                    
  real(8), save :: lay_height           (n_panel)                    
  real(8), save :: lay_x      (2, n_lay, n_panel)                    
  real(8), save :: lay_y      (2, n_lay, n_panel)                    
  real(8), save :: lay_center (2, n_lay, n_panel)                    
  real(8), save :: lay_velo   (2, n_lay, n_panel)                    
                                                                     
! *** Vortex element (blob) ***                                      
                                                                     
  integer, save :: nvor_b                                            
  integer, save :: blob_id              (n_vor_b)                    
                                                                     
  real(8), save :: vor_b             (2, n_vor_b)                    
  real(8), save :: vor_vb         (2, 3, n_vor_b)                    
  real(8), save :: cir_b                (n_vor_b)                    
  real(8), save :: cor_b                (n_vor_b)                    
                                                                     
! *** Vortex element (sheet) ***                                     
                                                                     
  integer, save :: nvor_s                                            
  integer, save :: sheet_id             (n_vor_s)                    
                                                                     
  real(8), save :: vor_s          (2, 2, n_vor_s)                    
  real(8), save :: vor_sc            (2, n_vor_s)                    
  real(8), save :: vor_vs         (2, 3, n_vor_s)                    
  real(8), save :: vor_sr               (n_vor_s)                    
  real(8), save :: cir_s                (n_vor_s)                    
  real(8), save :: cor_s                (n_vor_s)                    
                                                                     
! *** elapsed time check ***                                         
                                                                     
  integer, save :: hz                                                
  integer, save :: clock                                             
                                                                     
  real(8), save :: total_time                                        
  real(8), save :: time_matrix_velo                                  
  real(8), save :: time_matrix_pres                                  
  real(8), save :: time_get_velo                                     
  real(8), save :: time_core                                         
  real(8), save :: time_drift                                        
  real(8), save :: time_nascent                                      
  real(8), save :: time_trans                                        
                                                                     
! *** parallel computation ***                                       
                                                                     
  integer, save :: my_rank                                           
  integer, save :: num_procs                                         
  integer, save :: ierr                                              
  integer, save :: fquit                                             
                                                                     
! *** post processing ***                                            
                                                                     
  integer, save :: n_height                                          
  integer, save :: nx                                                
  integer, save :: ny                                                
  integer, save :: p_start                                           
  integer, save :: p_end                                             
  integer, save :: p_int                                             
                                                                     
  real(8), save :: cal_height                                        
  real(8), save :: cut_rd                                            
  real(8), save :: x_max                                             
  real(8), save :: x_min                                             
  real(8), save :: y_max                                             
  real(8), save :: y_min                                             
  real(8), save :: point    (2, n_field, n_field)                    
  real(8), save :: velo     (2, n_field, n_field)
             
                                                                     
 end module global                                                   
!====================================================================

