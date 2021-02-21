    
module MSTG_SEM_Lib
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    interface
        subroutine STG_compute_SEM_matrix_data_C( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33 &
        ) &
        bind(C, name='STG_compute_SEM_matrix_data')
            use, intrinsic :: ISO_C_BINDING
            real(c_float), intent(in), value ::  re_uu, re_vv, re_ww, re_uv, re_uw, re_vw
            type(c_ptr), value :: a11, a12, a13
            type(c_ptr), value :: a21, a22, a23
            type(c_ptr), value :: a31, a32, a33
        end subroutine
    end interface 
    
    interface 
        subroutine STG_compute_SEM_in_planes_lims_C( &
            x_min, x_max, y_min, y_max, z_min, z_max, &
	        x_e, y_e, z_e, num_eddies, &
	        u_e, v_e, w_e, &
	        x_min_in, x_max_in, &
	        y_min_in, y_max_in, &
	        z_min_in, z_max_in &
        ) &
        bind(C, name='STG_compute_SEM_in_planes_lims_fort')
            use, intrinsic :: ISO_C_BINDING
            real(c_float), intent(in), value ::  x_min, x_max, y_min, y_max, z_min, z_max
            type(c_ptr), value :: x_e, y_e, z_e 
            real(c_float), intent(in), value :: u_e, v_e, w_e
            integer(c_long), intent(in), value :: num_eddies
            type(c_ptr), value :: x_min_in, x_max_in, y_min_in, y_max_in, z_min_in, z_max_in
        end subroutine
    end interface
    
    interface
        subroutine STG_compute_SEM_init_eddies_params_C( &
            x_e, y_e, z_e, &
            eps_x, eps_y, eps_z, num_eddies, &
            x_min, x_max, &
            y_min, y_max, &
            z_min, z_max &
        ) &
        bind(C, name='STG_compute_SEM_init_eddies_params_fort')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), value :: x_e, y_e, z_e, eps_x, eps_y, eps_z
            integer(c_long), intent(in), value :: num_eddies
            real(c_float), intent(in), value :: x_min, x_max, y_min, y_max, z_min, z_max
        end subroutine
    end interface 
    
    interface 
        subroutine STG_compute_SEM_new_eddies_params_C( &
            x_e_cur, y_e_cur, z_e_cur, &
            eps_x_cur, eps_y_cur, eps_z_cur, num_eddies, &
            x_min, x_max, &
            y_min, y_max, &
            z_min, z_max, &
            u_e, v_e, w_e, ts, &
            x_min_in, x_max_in, &
            y_min_in, y_max_in, &
            z_min_in, z_max_in, &
            x_e_new, y_e_new, z_e_new, &
            eps_x_new, eps_y_new, eps_z_new &
        ) &
        bind(C, name='STG_compute_SEM_new_eddies_params_fort')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), value :: x_e_cur, y_e_cur, z_e_cur, eps_x_cur, eps_y_cur, eps_z_cur
            integer(c_long), value :: num_eddies
            real(c_float), intent(in), value :: x_min, x_max, y_min, y_max, z_min, z_max
            real(c_float), intent(in), value :: u_e, v_e, w_e, ts
            type(c_ptr), value :: x_min_in, x_max_in, y_min_in, y_max_in,z_min_in, z_max_in
            type(c_ptr), value :: x_e_new, y_e_new, z_e_new
            type(c_ptr), value :: eps_x_new, eps_y_new, eps_z_new
        end subroutine
    end interface
    
    interface 
        subroutine STG_compute_SEM_pulsations_C( &
            x_e, y_e, z_e, eps_x, eps_y, eps_z, &
            num_eddies, x, y, z, volume, &
            ls_ux, ls_uy, ls_uz, &
            ls_vx, ls_vy, ls_vz, &
            ls_wx, ls_wy, ls_wz, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            u, v, w &
        ) &
        bind(C, name='STG_compute_SEM_pulsations')
            use, intrinsic :: ISO_C_BINDING
            type(c_ptr), value :: x_e, y_e, z_e, eps_x, eps_y, eps_z
            integer(c_long), intent(in), value :: num_eddies
            real(c_float), intent(in), value :: x, y, z, volume
            real(c_float), intent(in), value :: ls_ux, ls_uy, ls_uz
            real(c_float), intent(in), value :: ls_vx, ls_vy, ls_vz
            real(c_float), intent(in), value :: ls_wx, ls_wy, ls_wz
            real(c_float), intent(in), value :: a11, a12, a13
            real(c_float), intent(in), value :: a21, a22, a23
            real(c_float), intent(in), value :: a31, a32, a33
            type(c_ptr), value :: u, v, w
        end subroutine
    end interface
    
    
    
    contains 
    
    
    subroutine STG_compute_SEM_matrix_data( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33 &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in) ::  re_uu, re_vv, re_ww, re_uv, re_uw, re_vw
        real :: a11, a12, a13
        real :: a21, a22, a23
        real :: a31, a32, a33
        
        call STG_compute_SEM_matrix_data_C( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            c_loc(a11), c_loc(a12), c_loc(a13), &
            c_loc(a21), c_loc(a22), c_loc(a23), &
            c_loc(a31), c_loc(a32), c_loc(a33) &
        )
    end subroutine
    
    
    subroutine STG_compute_SEM_in_planes_lims( &
        x_min, x_max, y_min, y_max, z_min, z_max, &
        x_e, y_e, z_e, num_eddies, &
        u_e, v_e, w_e, &
        x_min_in, x_max_in, &
        y_min_in, y_max_in, &
        z_min_in, z_max_in &
    )
        use, intrinsic :: ISO_C_BINDING
        real, intent(in) ::  x_min, x_max, y_min, y_max, z_min, z_max
        real :: x_e(:), y_e(:), z_e(:) 
        real, intent(in) :: u_e, v_e, w_e
        integer, intent(in) :: num_eddies
        real :: x_min_in(:), x_max_in(:), y_min_in(:), y_max_in(:), z_min_in(:), z_max_in(:)
        
        call STG_compute_SEM_in_planes_lims_C( &
            x_min, x_max, y_min, y_max, z_min, z_max, &
            c_loc(x_e(1)), c_loc(y_e(1)), c_loc(z_e(1)), num_eddies, &
            u_e, v_e, w_e, &
            c_loc(x_min_in(1)), c_loc(x_max_in(1)), &
            c_loc(y_min_in(1)), c_loc(y_max_in(1)), &
            c_loc(z_min_in(1)), c_loc(z_max_in(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_SEM_init_eddies_params( &
        x_e, y_e, z_e, &
        eps_x, eps_y, eps_z, num_eddies, &
        x_min, x_max, &
        y_min, y_max, &
        z_min, z_max &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: x_e(:), y_e(:), z_e(:), eps_x(:), eps_y(:), eps_z(:)
        integer, intent(in) :: num_eddies
        real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
        
        call STG_compute_SEM_init_eddies_params_C( &
            c_loc(x_e(1)), c_loc(y_e(1)), c_loc(z_e(1)), &
            c_loc(eps_x(1)), c_loc(eps_y(1)), c_loc(eps_z(1)), num_eddies, &
            x_min, x_max, &
            y_min, y_max, &
            z_min, z_max &
        )
    end subroutine
    
    
    
    subroutine STG_compute_SEM_new_eddies_params( &
        x_e_cur, y_e_cur, z_e_cur, &
        eps_x_cur, eps_y_cur, eps_z_cur, num_eddies, &
        x_min, x_max, &
        y_min, y_max, &
        z_min, z_max, &
        u_e, v_e, w_e, ts, &
        x_min_in, x_max_in, &
        y_min_in, y_max_in, &
        z_min_in, z_max_in, &
        x_e_new, y_e_new, z_e_new, &
        eps_x_new, eps_y_new, eps_z_new &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: x_e_cur(:), y_e_cur(:), z_e_cur(:), eps_x_cur(:), eps_y_cur(:), eps_z_cur(:)
        integer :: num_eddies
        real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
        real, intent(in) :: u_e, v_e, w_e, ts
        real :: x_min_in(:), x_max_in(:), y_min_in(:), y_max_in(:), z_min_in(:), z_max_in(:)
        real :: x_e_new(:), y_e_new(:), z_e_new(:)
        real :: eps_x_new(:), eps_y_new(:), eps_z_new(:)
        
        call STG_compute_SEM_new_eddies_params_C( &
            c_loc(x_e_cur(1)), c_loc(y_e_cur(1)), c_loc(z_e_cur(1)), &
            c_loc(eps_x_cur(1)), c_loc(eps_y_cur(1)), c_loc(eps_z_cur(1)), num_eddies, &
            x_min, x_max, &
            y_min, y_max, &
            z_min, z_max, &
            u_e, v_e, w_e, ts, &
            c_loc(x_min_in(1)), c_loc(x_max_in(1)), &
            c_loc(y_min_in(1)), c_loc(y_max_in(1)), &
            c_loc(z_min_in(1)), c_loc(z_max_in(1)), &
            c_loc(x_e_new(1)), c_loc(y_e_new(1)), c_loc(z_e_new(1)), &
            c_loc(eps_x_new(1)), c_loc(eps_y_new(1)), c_loc(eps_z_new(1)) &
        )
    end subroutine
    
    
    subroutine STG_compute_SEM_pulsations( &
        x_e, y_e, z_e, eps_x, eps_y, eps_z, &
        num_eddies, x, y, z, volume, &
        ls_ux, ls_uy, ls_uz, &
        ls_vx, ls_vy, ls_vz, &
        ls_wx, ls_wy, ls_wz, &
        a11, a12, a13, &
        a21, a22, a23, &
        a31, a32, a33, &
        u, v, w &
    )
        use, intrinsic :: ISO_C_BINDING
        real :: x_e(:), y_e(:), z_e(:), eps_x(:), eps_y(:), eps_z(:)
        integer, intent(in) :: num_eddies
        real, intent(in) :: x, y, z, volume
        real, intent(in) :: ls_ux, ls_uy, ls_uz
        real, intent(in) :: ls_vx, ls_vy, ls_vz
        real, intent(in) :: ls_wx, ls_wy, ls_wz
        real, intent(in) :: a11, a12, a13
        real, intent(in) :: a21, a22, a23
        real, intent(in) :: a31, a32, a33
        real, intent(out) :: u, v, w
        
        call STG_compute_SEM_pulsations_C( &
            c_loc(x_e(1)), c_loc(y_e(1)), c_loc(z_e(1)), c_loc(eps_x(1)), c_loc(eps_y(1)), c_loc(eps_z(1)), &
            num_eddies, x, y, z, volume, &
            ls_ux, ls_uy, ls_uz, &
            ls_vx, ls_vy, ls_vz, &
            ls_wx, ls_wy, ls_wz, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            c_loc(u), c_loc(v), c_loc(w) &
        )
    end subroutine
    
end module MSTG_SEM_Lib