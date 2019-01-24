    
module STG_SMIRNOV
    
    use, intrinsic :: ISO_C_BINDING
    implicit none

    
    interface 
        subroutine STG_compute_Smirnov_matrix_data_homo( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33 &
        ) & 
        bind(C, name='STG_compute_Smirnov_matrix_data_homo')
            use ISO_C_BINDING
            real(c_float), intent(in), value :: re_uu, re_vv, re_ww
            real(c_float), intent(in), value :: re_uv, re_uw, re_vw
            type(c_ptr), value :: c1, c2, c3
            type(c_ptr), value :: a11, a12, a13
            type(c_ptr), value :: a21, a22, a23
            type(c_ptr), value :: a31, a32, a33
        end subroutine STG_compute_Smirnov_matrix_data_homo
    end interface 
    
    interface 
        subroutine STG_compute_Smirnov_random_data( &
            num_modes, omega, &
            k1, k2, k3, &
            zeta1, zeta2, zeta3, &
            xi1, xi2, xi3, &
            p1, p2, p3, &
            q1, q2, q3 &
        ) & 
        bind(C, name='STG_compute_Smirnov_random_data')
            use ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: k1, k2, k3, omega
            type(c_ptr), value :: zeta1, zeta2, zeta3, xi1, xi2, xi3
            type(c_ptr), value :: p1, p2, p3, q1, q2, q3
        end subroutine STG_compute_Smirnov_random_data
    end interface 
    
end module STG_SMIRNOV