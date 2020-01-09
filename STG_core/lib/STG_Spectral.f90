    
module STG_SPECTRAL
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    
    interface 
        subroutine STG_compute_Spectral_random_angles_and_phase_C( &
            num_modes, phi, psi, &
            alpha, theta, omega &
        ) &
        bind(C, name='STG_compute_Spectral_random_angles_and_phase')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: phi, psi, alpha, theta, omega
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Spectral_modes_params_C( &
            num_modes, k_arr, phi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
        ) &
        bind(C, name='STG_compute_Spectral_modes_params')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: k_arr, phi, alpha, theta
            type(c_ptr), value :: k1, k2, k3, sigma1, sigma2, sigma3
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Spectral_matrix_data_C( &
            re_uu, re_vv, re_ww, &
	        re_uv, re_uw, re_vw, &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33  &
        ) &
        bind(C, name='STG_compute_Spectral_matrix_data')
            use, intrinsic :: ISO_C_BINDING
            real(c_float), intent(in), value :: re_uu, re_vv, re_ww
            real(c_float), intent(in), value :: re_uv, re_uw, re_vw
            type(c_ptr), value :: c1, c2, c3
            type(c_ptr), value :: a11, a12, a13
            type(c_ptr), value :: a21, a22, a23
            type(c_ptr), value :: a31, a32, a33
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Spectral_pulsations_C( &
            k1, k2, k3, &
	        sigma1, sigma2, sigma3, psi, omega, u_abs, &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33, &
	        x, y, z, &
	        num_modes, time_scale, time, &
	        u_p, v_p, w_p &
        ) &
        bind(C, name='STG_compute_Spectral_pulsations')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), value :: k1, k2, k3, sigma1, sigma2, sigma3, psi, omega, u_abs
            real(c_float), intent(in), value :: c1, c2, c3, a11, a12, a13
            real(c_float), intent(in), value :: a21, a22, a23, a31, a32, a33
            real(c_float), intent(in), value :: x, y, z, time_scale, time
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: u_p, v_p, w_p 
        end subroutine
    end interface
    
    
    
    
    contains
    
    
    subroutine STG_compute_Spectral_random_angles_and_phase( &
        num_modes, phi, psi, &
        alpha, theta, omega &
    )
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: num_modes
        real ::phi(:), psi(:), alpha(:), theta(:), omega(:)
        
        call STG_compute_Spectral_random_angles_and_phase_C( &
            num_modes, c_loc(phi(1)), c_loc(psi(1)), &
            c_loc(alpha(1)), c_loc(theta(1)), c_loc(omega(1)) &
        )
    end subroutine
    
    
    
    subroutine STG_compute_Spectral_modes_params( &
        num_modes, k_arr, phi, &
        alpha, theta, k1, k2, k3, &
        sigma1, sigma2, sigma3 &
    )
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: num_modes
        real :: k_arr(:), phi(:), alpha(:), theta(:)
        real :: k1(:), k2(:), k3(:), sigma1(:), sigma2(:), sigma3(:)
        
        call STG_compute_Spectral_modes_params_C( &
            num_modes, c_loc(k_arr(1)), c_loc(phi(1)), &
            c_loc(alpha(1)), c_loc(theta(1)), & 
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(sigma1(1)), c_loc(sigma2(1)), c_loc(sigma3(1)) &
        )
    end subroutine
    
    
    
    subroutine STG_compute_Spectral_matrix_data( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        c1, c2, c3, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33  &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in) :: re_uu, re_vv, re_ww 
        real, intent(in) :: re_uv, re_uw, re_vw
        real :: c1, c2, c3, a11, a12, a13
        real :: a21, a22, a23, a31, a32, a33
        
        call STG_compute_Spectral_matrix_data_C( &
            re_uu, re_vv, re_ww, &
	        re_uv, re_uw, re_vw, &
	        c_loc(c1), c_loc(c2), c_loc(c3), &
	        c_loc(a11), c_loc(a12), c_loc(a13), &
	        c_loc(a21), c_loc(a22), c_loc(a23), &
	        c_loc(a31), c_loc(a32), c_loc(a33)  &
        )
    end subroutine
    
    
    
    subroutine STG_compute_Spectral_pulsations( &
        k1, k2, k3, &
        sigma1, sigma2, sigma3, psi, omega, u_abs, &
        c1, c2, c3, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33, &
        x, y, z, &
        num_modes, time_scale, time, &
        u_p, v_p, w_p &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: k1(:), k2(:), k3(:)
        real :: sigma1(:), sigma2(:), sigma3(:), psi(:), omega(:), u_abs(:)
        real, intent(in) :: c1, c2, c3, a11, a12, a13
        real, intent(in) :: a21, a22, a23, a31, a32, a33
        real, intent(in) :: x, y, z, time_scale, time
        integer, intent(in) :: num_modes
        real, intent(in) :: u_p, v_p, w_p
        
        call STG_compute_Spectral_pulsations_C( &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
	        c_loc(sigma1(1)), c_loc(sigma2(1)), c_loc(sigma3(1)), & 
            c_loc(psi(1)), c_loc(omega(1)), c_loc(u_abs(1)), &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33, &
	        x, y, z, &
	        num_modes, time_scale, time, &
	        c_loc(u_p), c_loc(v_p), c_loc(w_p) &
        )
    end subroutine
    
   
end module STG_SPECTRAL