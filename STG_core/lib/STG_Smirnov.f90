    
module STG_SMIRNOV
    
    use, intrinsic :: ISO_C_BINDING
    implicit none

    interface 
        subroutine STG_compute_Smirnov_matrix_data_homo_C( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33 &
        ) & 
        bind(C, name='STG_compute_Smirnov_matrix_data_homo')
            use, intrinsic :: ISO_C_BINDING
            real(c_float), intent(in), value :: re_uu, re_vv, re_ww
            real(c_float), intent(in), value :: re_uv, re_uw, re_vw
            type(c_ptr), value :: c1, c2, c3
            type(c_ptr), value :: a11, a12, a13
            type(c_ptr), value :: a21, a22, a23
            type(c_ptr), value :: a31, a32, a33
        end subroutine STG_compute_Smirnov_matrix_data_homo_C
    end interface 
    
    interface 
        subroutine STG_compute_Smirnov_random_data_C( &
            num_modes, omega, &
            k1, k2, k3, &
            zeta1, zeta2, zeta3, &
            xi1, xi2, xi3, &
            p1, p2, p3, &
            q1, q2, q3 &
        ) & 
        bind(C, name='STG_compute_Smirnov_random_data')
            use, intrinsic :: ISO_C_BINDING
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: k1, k2, k3, omega
            type(c_ptr), value :: zeta1, zeta2, zeta3, xi1, xi2, xi3
            type(c_ptr), value :: p1, p2, p3, q1, q2, q3
        end subroutine STG_compute_Smirnov_random_data_C
    end interface
    
    interface 
        subroutine STG_compute_Smirnov_pulsations_C( &
            k1, k2, k3, &
            p1, p2, p3, &
            q1, q2, q3, &
            omega, &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            x, y, z, &
            length_scale, time_scale, &
            num_modes, time, &
            u, v, w &
        ) & 
        bind(C, name='STG_compute_Smirnov_pulsations')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), intent(in), value :: k1, k2, k3, omega
            type(c_ptr), intent(in), value :: p1, p2, p3, q1, q2, q3
            real(c_float), intent(in), value :: c1, c2, c3
            real(c_float), intent(in), value :: a11, a12, a13
            real(c_float), intent(in), value :: a21, a22, a23
            real(c_float), intent(in), value :: a31, a32, a33
            integer(c_long), intent(in), value :: num_modes
            real(c_float), intent(in), value :: x, y, z, length_scale, time_scale, time
            type(c_ptr), value :: u, v, w
        end subroutine STG_compute_Smirnov_pulsations_C
    end interface
    
    contains
    
    subroutine STG_compute_Smirnov_matrix_data_homo( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        c1, c2, c3, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33 &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in), value :: re_uu, re_vv, re_ww
        real, intent(in), value :: re_uv, re_uw, re_vw
        real :: c1, c2, c3
        real :: a11, a12, a13
        real :: a21, a22, a23
        real :: a31, a32, a33
        
        call STG_compute_Smirnov_matrix_data_homo_C( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            c_loc(c1), c_loc(c2), c_loc(c3), &
            c_loc(a11), c_loc(a12), c_loc(a13), &
            c_loc(a21), c_loc(a22), c_loc(a23), &
            c_loc(a31), c_loc(a32), c_loc(a33) &
        )
    end subroutine STG_compute_Smirnov_matrix_data_homo
        
    subroutine STG_compute_Smirnov_random_data( &
        num_modes, omega, &
        k1, k2, k3, &
        zeta1, zeta2, zeta3, &
        xi1, xi2, xi3, &
        p1, p2, p3, &
        q1, q2, q3 &
    )
        use, intrinsic :: ISO_C_BINDING
        integer(8), intent(in) :: num_modes
        real :: k1(:), k2(:), k3(:), omega(:)
        real :: zeta1(:), zeta2(:), zeta3(:), xi1(:), xi2(:), xi3(:)
        real :: p1(:), p2(:), p3(:), q1(:), q2(:), q3(:)
        
        call STG_compute_Smirnov_random_data_C( &
            num_modes, c_loc(omega(1)), &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(zeta1(1)), c_loc(zeta2(1)), c_loc(zeta3(1)), &
            c_loc(xi1(1)), c_loc(xi2(1)), c_loc(xi3(1)), &
            c_loc(p1(1)), c_loc(p2(1)), c_loc(p3(1)), &
            c_loc(q1(1)), c_loc(q2(1)), c_loc(q3(1)) &
        )
    end subroutine STG_compute_Smirnov_random_data
    
    subroutine STG_compute_Smirnov_pulsations( &
        k1, k2, k3, &
        p1, p2, p3, &
        q1, q2, q3, &
        omega, &
        c1, c2, c3, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33, &
        x, y, z, &
        length_scale, time_scale, &
        num_modes, time, &
        u, v, w &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in) :: k1(:), k2(:), k3(:), omega(:)
        real, intent(in) :: p1(:), p2(:), p3(:), q1(:), q2(:), q3(:)
        real, intent(in) :: c1, c2, c3
        real, intent(in) :: a11, a12, a13
        real, intent(in) :: a21, a22, a23
        real, intent(in) :: a31, a32, a33
        integer(8), intent(in) :: num_modes
        real, intent(in) :: x, y, z, length_scale, time_scale, time
        real, intent(out) :: u, v, w
        
        call STG_compute_Smirnov_pulsations_C( &
            c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
            c_loc(p1(1)), c_loc(p2(1)), c_loc(p3(1)), &
            c_loc(q1(1)), c_loc(q2(1)), c_loc(q3(1)), &
            c_loc(omega(1)), &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            x, y, z, &
            length_scale, time_scale, &
            num_modes, time, &
            c_loc(u), c_loc(v), c_loc(w) &
        )
    end subroutine STG_compute_Smirnov_pulsations
        
    
end module STG_SMIRNOV