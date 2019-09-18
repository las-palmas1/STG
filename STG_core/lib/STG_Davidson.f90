    
module STG_DAVIDSON
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    
    interface 
        subroutine STG_compute_Davidson_random_angles_and_phase_C( &
            num_modes, phi, psi, &
            alpha, theta &
        ) &
        bind(C, name='STG_compute_Davidson_random_angles_and_phase')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: phi, psi, alpha, theta
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Davidson_modes_params_C( &
            num_modes, k_arr, phi, psi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
        ) &
        bind(C, name='STG_compute_Davidson_modes_params')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: k_arr, phi, psi, alpha, theta
            type(c_ptr), value :: k1, k2, k3, sigma1, sigma2, sigma3
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Davidson_random_data_C( &
            num_modes, k_arr, phi, psi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
        ) &
        bind(C, name='STG_compute_Davidson_random_data')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: k_arr, phi, psi, alpha, theta
            type(c_ptr), value :: k1, k2, k3, sigma1, sigma2, sigma3
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Davidson_spectrum_C( &
            delta_min, num_modes, re_uu, re_vv, re_ww, &
            ls_i, dissip_rate, visc, energy, k_arr, u_abs &
        ) &
        bind(C, name='STG_compute_Davidson_spectrum')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            real(c_float), intent(in), value :: delta_min, re_uu, re_vv, re_ww
            real(c_float), intent(in), value :: ls_i, dissip_rate, visc
            type(c_ptr), value :: energy, k_arr, u_abs
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Davidson_matrix_data_C( &
            re_uu, re_vv, re_ww, &
	        re_uv, re_uw, re_vw, &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33  &
        ) &
        bind(C, name='STG_compute_Davidson_matrix_data')
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
        subroutine STG_compute_Davidson_auto_coef_C( &
            ts, ts_i, a, b &
        ) &
        bind(C, name='STG_compute_Davidson_auto_coef')
            use, intrinsic :: ISO_C_BINDING
            real(c_float), intent(in), value :: ts, ts_i
            type(c_ptr), value :: a, b
        end subroutine
    end interface
    
    
    interface 
        subroutine STG_compute_Davidson_pulsations_C( &
            k1, k2, k3, &
	        sigma1, sigma2, sigma3, psi, u_abs, &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33, &
	        x, y, z, &
	        a, b, num_modes, &
	        u_p_prev, v_p_prev, w_p_prev, &
	        u_p, v_p, w_p &
        ) &
        bind(C, name='STG_compute_Davidson_pulsations')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), value :: k1, k2, k3, sigma1, sigma2, sigma3, psi, u_abs
            real(c_float), intent(in), value :: c1, c2, c3, a11, a12, a13
            real(c_float), intent(in), value :: a21, a22, a23, a31, a32, a33
            real(c_float), intent(in), value :: x, y, z, a, b
            integer(c_long), intent(in), value :: num_modes
            real(c_float), intent(in), value :: u_p_prev, v_p_prev, w_p_prev
            type(c_ptr), value :: u_p, v_p, w_p 
        end subroutine
    end interface
    
    
    contains
    
    
    subroutine STG_compute_Davidson_random_angles_and_phase( &
            num_modes, phi, psi, &
            alpha, theta &
    )
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: num_modes
        real ::phi(:), psi(:), alpha(:), theta(:)
            
        call STG_compute_Davidson_random_angles_and_phase_C( &
            num_modes, c_loc(phi(1)), c_loc(psi(1)), &
            c_loc(alpha(1)), c_loc(theta(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_Davidson_modes_params( &
            num_modes, k_arr, phi, psi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
    )
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: num_modes
        real :: k_arr(:), phi(:), psi(:), alpha(:), theta(:)
        real :: k1(:), k2(:), k3(:), sigma1(:), sigma2(:), sigma3(:)
            
        call STG_compute_Davidson_modes_params_C( &
            num_modes, c_loc(k_arr(1)), c_loc(phi(1)), c_loc(psi(1)), &
            c_loc(alpha(1)), c_loc(theta(1)), &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(sigma1(1)), c_loc(sigma2(1)), c_loc(sigma3(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_Davidson_random_data( &
            num_modes, k_arr, phi, psi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
    )
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: num_modes
        real :: k_arr(:), phi(:), psi(:), alpha(:), theta(:)
        real :: k1(:), k2(:), k3(:), sigma1(:), sigma2(:), sigma3(:)
            
        call STG_compute_Davidson_random_data_C( &
            num_modes, c_loc(k_arr(1)), c_loc(phi(1)), c_loc(psi(1)), &
            c_loc(alpha(1)), c_loc(theta(1)), &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(sigma1(1)), c_loc(sigma2(1)), c_loc(sigma3(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_Davidson_spectrum( &
        delta_min, num_modes, re_uu, re_vv, re_ww, &
        ls_i, dissip_rate, visc, energy, k_arr, u_abs &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in) :: delta_min
        integer, intent(in) :: num_modes
        real, intent(in) :: re_uu, re_vv, re_ww 
        real, intent(in) :: ls_i, dissip_rate, visc
        real :: energy(:), k_arr(:), u_abs(:) 
        
        call STG_compute_Davidson_spectrum_C( &
            delta_min, num_modes, re_uu, re_vv, re_ww, &
            ls_i, dissip_rate, visc, & 
            c_loc(energy(1)), c_loc(k_arr(1)), c_loc(u_abs(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_Davidson_matrix_data( &
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
        
        call STG_compute_Davidson_matrix_data_C( &
            re_uu, re_vv, re_ww, &
	        re_uv, re_uw, re_vw, &
	        c_loc(c1), c_loc(c2), c_loc(c3), &
	        c_loc(a11), c_loc(a12), c_loc(a13), &
	        c_loc(a21), c_loc(a22), c_loc(a23), &
	        c_loc(a31), c_loc(a32), c_loc(a33)  &
        )
    end subroutine
    
    
    subroutine STG_compute_Davidson_auto_coef( &
        ts, ts_i, a, b &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: ts, ts_i, a, b
        
        call STG_compute_Davidson_auto_coef_C( &
            ts, ts_i, c_loc(a), c_loc(b) &
        )
    end subroutine 
    
    
    subroutine STG_compute_Davidson_pulsations( &
        k1, k2, k3, &
        sigma1, sigma2, sigma3, psi, u_abs, &
        c1, c2, c3, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33, &
        x, y, z, &
        a, b, num_modes, &
        u_p_prev, v_p_prev, w_p_prev, &
        u_p, v_p, w_p &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: k1(:), k2(:), k3(:)
        real :: sigma1(:), sigma2(:), sigma3(:), psi(:), u_abs(:)
        real, intent(in) :: c1, c2, c3, a11, a12, a13
        real, intent(in) :: a21, a22, a23, a31, a32, a33
        real, intent(in) :: x, y, z, a, b
        integer, intent(in) :: num_modes
        real, intent(in) :: u_p_prev, v_p_prev, w_p_prev, u_p, v_p, w_p
        
        call STG_compute_Davidson_pulsations_C( &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(sigma1(1)), c_loc(sigma2(1)), c_loc(sigma3(1)), c_loc(psi(1)), c_loc(u_abs(1)), &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            x, y, z, &
            a, b, num_modes, &
            u_p_prev, v_p_prev, w_p_prev, &
            c_loc(u_p), c_loc(v_p), c_loc(w_p) &
        )
    end subroutine
    
end module STG_DAVIDSON