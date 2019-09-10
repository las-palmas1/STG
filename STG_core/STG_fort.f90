!  STG_fort.f90 
!
!  FUNCTIONS:
!  STG_fort - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: STG_fort
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
module Test

contains
    subroutine test_Smirnov( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        num_modes, ls_i &
    )
        use STG_COMMON
        use STG_SMIRNOV
        implicit none
        
        real :: re_uu, re_vv, re_ww
        real :: re_uv, re_uw, re_vw
        integer(8) :: num_modes
        
        integer :: num_rand_vars, num_mat_vars
        real :: time, x, y, z, length_scale, time_scale, ls_i
        
        real, allocatable :: rand_data(:, :)
        real, allocatable :: mat_data(:)
        real, allocatable :: k1(:), k2(:), k3(:), zeta1(:), zeta2(:), zeta3(:), xi1(:), xi2(:), xi3(:)
        real, allocatable :: p1(:), p2(:), p3(:), q1(:), q2(:), q3(:), omega(:)
        real :: c1, c2, c3, a11, a12, a13
        real :: a21, a22, a23, a31, a32, a33
        real :: u, v, w
        
        num_rand_vars = 16
        num_mat_vars = 12
        time = 0
        x = 1
        y = 1
        z = 1
        length_scale = ls_i
        time_scale = 1
        
        allocate(rand_data(num_modes, num_rand_vars))
        allocate(k1(num_modes))
        allocate(k2(num_modes))
        allocate(k3(num_modes))
        allocate(zeta1(num_modes))
        allocate(zeta2(num_modes))
        allocate(zeta3(num_modes))
        allocate(xi1(num_modes))
        allocate(xi2(num_modes))
        allocate(xi3(num_modes))
        allocate(p1(num_modes))
        allocate(p2(num_modes))
        allocate(p3(num_modes))
        allocate(q1(num_modes))
        allocate(q2(num_modes))
        allocate(q3(num_modes))
        allocate(omega(num_modes))
        
        allocate(mat_data(num_mat_vars))
        
        call STG_init_rand()
        call STG_compute_Smirnov_random_data( &
            num_modes, omega, &
            k1, k2, k3, &
            zeta1, zeta2, zeta3, &
            xi1, xi2, xi3, &
            p1, p2, p3, &
            q1, q2, q3 &
        )
        call STG_compute_Smirnov_matrix_data_homo( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33 &
        )
            
        rand_data(:, 1) = k1
        rand_data(:, 2) = k2
        rand_data(:, 3) = k3
        rand_data(:, 4) = zeta1
        rand_data(:, 5) = zeta2
        rand_data(:, 6) = zeta3
        rand_data(:, 7) = xi1
        rand_data(:, 8) = xi2
        rand_data(:, 9) = xi3
        rand_data(:, 10) = p1
        rand_data(:, 11) = p2
        rand_data(:, 12) = p3
        rand_data(:, 13) = q1
        rand_data(:, 14) = q2
        rand_data(:, 15) = q3
        rand_data(:, 16) = omega
        
        mat_data(1) = c1
        mat_data(2) = c2
        mat_data(3) = c3
        mat_data(4) = a11
        mat_data(5) = a12
        mat_data(6) = a13
        mat_data(7) = a21
        mat_data(8) = a22
        mat_data(9) = a23
        mat_data(10) = a31
        mat_data(11) = a32
        mat_data(12) = a33
        
        call STG_compute_Smirnov_pulsations( &
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
        
        print *, '   Test Smirnov'
        print *, 'Reynolds stresses'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uu, re_uv, re_uw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uv, re_vv, re_vw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uw, re_vw, re_ww
        print '(1x, a4, f8.3, 2x, a4, f8.3, 2x, a4, f8.3)', 'u = ', u, 'v = ', v, 'w = ', w
        print *, ''
        
        deallocate(rand_data)
        deallocate(k1, k2, k3, zeta1, zeta2, zeta3)
        deallocate(xi1, xi2, xi3, p1, p2, p3, q1, q2, q3, omega)
    end subroutine test_Smirnov
    
    
    
    subroutine test_Davidson( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        num_modes, ls_i &
    )
        use STG_COMMON
        use STG_DAVIDSON
        real :: re_uu, re_vv, re_ww, re_uv, re_uw, re_vw
        integer :: num_modes
        real:: ls_i
        
        real :: dissip_rate, visc, delta_min
        real :: c1, c2, c3, a11, a12, a13
        real :: a21, a22, a23, a31, a32, a33
        real :: ts_i, ts, a, b, x, y, z
        real :: u, v, w, u_prev, v_prev, w_prev
        real, allocatable :: energy(:), k_arr(:), u_abs(:)
        real, allocatable :: phi(:), psi(:), alpha(:), theta(:)
        real, allocatable :: k1(:), k2(:), k3(:), sigma1(:), sigma2(:), sigma3(:)
        
        call STG_init_rand()
        call STG_compute_Davidson_matrix_data( &
            re_uu, re_vv, re_ww, &
	        re_uv, re_uw, re_vw, &
	        c1, c2, c3, &
	        a11, a12, a13, &
	        a21, a22, a23, &
	        a31, a32, a33  &
        )
        
        visc = 1.5e-5
        dissip_rate = 0.09**0.75 * (0.5 * (re_uu + re_vv + re_ww))**1.5 / ls_i
        delta_min = 0.001
        ts_i = 0.001
        ts = 0.01
        
        allocate(energy(num_modes))
        allocate(k_arr(num_modes))
        allocate(u_abs(num_modes))
        
        call STG_compute_Davidson_spectrum( &
            delta_min, num_modes, re_uu, re_vv, re_ww, &
            ls_i, dissip_rate, visc, energy, k_arr, u_abs &
        )
        
        allocate(phi(num_modes), psi(num_modes), alpha(num_modes), theta(num_modes))
        allocate(k1(num_modes), k2(num_modes), k3(num_modes))
        allocate(sigma1(num_modes), sigma2(num_modes), sigma3(num_modes))
        
        call STG_compute_Davidson_random_data( &
            num_modes, k_arr, phi, psi, &
            alpha, theta, k1, k2, k3, &
            sigma1, sigma2, sigma3 &
        )
            
        call STG_compute_Davidson_auto_coef(ts, ts_i, a, b)
        print *, 'a b'
        print '(1x, f11.8, 2x, f11.8)', a, b
        a = 0.
        b = 1.
        u_prev = 0.
        v_prev = 0.
        w_prev = 0.
        x = 1.
        y = 1.
        z = 1.
        
        call STG_compute_Davidson_pulsations( &
            k1, k2, k3, &
            sigma1, sigma2, sigma3, psi, u_abs, &
            c1, c2, c3, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33, &
            x, y, z, &
            a, b, num_modes, &
            u_prev, v_prev, w_prev, &
            u, v, w &
        )
            
        print *, '   Test Davidson'
        print *, 'Reynolds stresses'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uu, re_uv, re_uw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uv, re_vv, re_vw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uw, re_vw, re_ww
        print *, 'Eig values'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', c1, c2, c3
        print *, ' U_abs'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2, 2x, f6.2)', u_abs(1: 4)
        print *, 'Rand data (theta)'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2, 2x, f6.2)', theta(1: 4) 
        print *, 'Velocities'
        print '(1x, a4, f8.3, 2x, a4, f8.3, 2x, a4, f8.3)', 'u = ', u, 'v = ', v, 'w = ', w
        print *, ''
        
        deallocate(u_abs, k_arr, energy)
        deallocate(phi, psi, alpha, theta)
        deallocate(k1, k2, k3)
        deallocate(sigma1, sigma2, sigma3)
    end subroutine 
    
    
    subroutine test_SEM( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        num_eddies, & 
        u_e, v_e, w_e, &
        ls_ux, ls_uy, ls_uz, &
        ls_vx, ls_vy, ls_vz, &
        ls_wx, ls_wy, ls_wz &
    )
        use STG_COMMON
        use STG_SEM
        
        real :: re_uu, re_vv, re_ww, re_uv, re_uw, re_vw
        integer :: num_eddies
        real :: u_e, v_e, w_e
        real:: ls_ux, ls_uy, ls_uz
        real:: ls_vx, ls_vy, ls_vz
        real:: ls_wx, ls_wy, ls_wz
        
        real :: a11, a12, a13
        real :: a21, a22, a23
        real :: a31, a32, a33
        real :: x_min, x_max, y_min, y_max, z_min, z_max, volume, x, y, z
        real, allocatable :: x_e(:), y_e(:), z_e(:), eps_x(:), eps_y(:), eps_z(:)
        real, allocatable :: x_e_new(:), y_e_new(:), z_e_new(:), eps_x_new(:), eps_y_new(:), eps_z_new(:)
        real, allocatable :: x_min_in(:), x_max_in(:), y_min_in(:), y_max_in(:), z_min_in(:), z_max_in(:)
        real :: ts
        real :: u, v, w
        
        ts = 0.5
        x_min = 0.
        y_min = 0.
        z_min = 0.
        x_max = 1.
        y_max = 1.
        z_max = 1.
        volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
        
        allocate(x_e(num_eddies), y_e(num_eddies), z_e(num_eddies)) 
        allocate(eps_x(num_eddies), eps_y(num_eddies), eps_z(num_eddies)) 
        allocate(x_e_new(num_eddies), y_e_new(num_eddies), z_e_new(num_eddies))
        allocate(eps_x_new(num_eddies), eps_y_new(num_eddies), eps_z_new(num_eddies))  
        allocate(x_min_in(num_eddies), x_max_in(num_eddies))
        allocate(y_min_in(num_eddies), y_max_in(num_eddies))
        allocate(z_min_in(num_eddies), z_max_in(num_eddies))
        
        call STG_init_rand() 
        
        call STG_compute_SEM_matrix_data( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            a11, a12, a13, &
            a21, a22, a23, &
            a31, a32, a33 &
        )
            
        call STG_compute_SEM_init_eddies_params( &
            x_e, y_e, z_e, &
            eps_x, eps_y, eps_z, num_eddies, &
            x_min, x_max, &
            y_min, y_max, &
            z_min, z_max &
        )
            
        call STG_compute_SEM_in_planes_lims( &
            x_min, x_max, y_min, y_max, z_min, z_max, &
            x_e, y_e, z_e, num_eddies, &
            u_e, v_e, w_e, &
            x_min_in, x_max_in, &
            y_min_in, y_max_in, &
            z_min_in, z_max_in &
        )
            
        call STG_compute_SEM_new_eddies_params( &
            x_e, y_e, z_e, &
            eps_x, eps_y, eps_z, num_eddies, &
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
            
        x = 0.5
        y = 0.5
        z = 0.5
            
        call STG_compute_SEM_pulsations( &
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
        
        print *, '   Test SEM '
        print *, 'Reynolds stresses'
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uu, re_uv, re_uw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uv, re_vv, re_vw
        print '(1x, f6.2, 2x, f6.2, 2x, f6.2)', re_uw, re_vw, re_ww
        print *, 'Matrix data'
        print '(1x, a4, f6.2, a4, f6.2, a4, f6.2)', ' a11= ', a11, ' a12= ', a12, ' a13= ', a13
        print '(1x, a4, f6.2, a4, f6.2, a4, f6.2)', ' a21= ', a21, ' a22= ', a22, ' a23= ', a23
        print '(1x, a4, f6.2, a4, f6.2, a4, f6.2)', ' a31= ', a31, ' a32= ', a32, ' a33= ', a33
        print *, 'Init eddies params'
        print '(1x, a8, f6.3, a8, f6.3, a8, f6.3)', ' x_e= ', x_e(1), ' y_e= ', y_e(1), ' z_e= ', z_e(1)
        print '(1x, a8, f6.3, a8, f6.3, a8, f6.3)', ' eps_x= ', eps_x(1), ' eps_y= ', eps_y(1), ' eps_z= ', eps_z(1)
        print *, 'In planes lims'
        print '(1x, a12, f6.2, a12, f6.2, a12, f6.2)', ' x_min_in= ', x_min_in(1), ' y_min_in= ', y_min_in(1), ' z_min_in= ', z_min_in(1)
        print '(1x, a12, f6.2, a12, f6.2, a12, f6.2)', ' x_max_in= ', x_max_in(1), ' y_max_in= ', y_max_in(1), ' z_max_in= ', z_max_in(1)
        print *, 'New eddies params'
        print '(1x, a8, f6.3, a8, f6.3, a8, f6.3)', ' x_e= ', x_e_new(1), ' y_e= ', y_e_new(1), ' z_e= ', z_e_new(1)
        print '(1x, a8, f6.3, a8, f6.3, a8, f6.3)', ' eps_x= ', eps_x_new(1), ' eps_y= ', eps_y_new(1), ' eps_z= ', eps_z_new(1)
        print *, 'Velocities'
        print '(1x, a4, f8.3, 2x, a4, f8.3, 2x, a4, f8.3)', 'u = ', u, 'v = ', v, 'w = ', w
        print *, ''
        
        if (ANY(ISNAN(x_e))) print *, 'x_e NAN'
        if (ANY(ISNAN(y_e))) print *, 'y_e NAN'
        if (ANY(ISNAN(z_e))) print *, 'z_e NAN'
        if (ANY(ISNAN(eps_x))) print *, 'eps_x NAN'
        if (ANY(ISNAN(eps_y))) print *, 'eps_y NAN'
        if (ANY(IsNAN(eps_z))) print *, 'eps_z NAN'
        
        deallocate(x_e, y_e, z_e)
        deallocate(eps_x, eps_y, eps_z)
        deallocate(x_e_new, y_e_new, z_e_new)
        deallocate(eps_x_new, eps_y_new, eps_z_new)  
        deallocate(x_min_in, x_max_in)
        deallocate(y_min_in, y_max_in)
        deallocate(z_min_in, z_max_in)
    
    end subroutine
    
    
    subroutine Read_array(turb_data_dir, dir, fname, arr)
        character(len=100) :: turb_data_dir, full_fname
        character(*) :: fname, dir
        real :: arr(:)
     
        full_fname = TRIM(turb_data_dir) // "\" // TRIM(dir) // "\" // TRIM(fname)
        open(21, file = full_fname) 
        read(21, *) arr
        close(21)
    end subroutine 

    
    subroutine Read_re_stresses(turb_data_dir, y, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw)
        character(len=100) :: turb_data_dir
        real :: y(:), re_uu(:), re_vv(:), re_ww(:), re_uv(:), re_uw(:), re_vw(:)
    
        call Read_array(turb_data_dir, "re", "uu_av", re_uu)
        call Read_array(turb_data_dir, 're', 'vv_av', re_vv)
        call Read_array(turb_data_dir, 're', 'ww_av', re_ww)
        call Read_array(turb_data_dir, 're', 'uv_av', re_uv)
        call Read_array(turb_data_dir, 're', 'uw_av', re_uw)
        call Read_array(turb_data_dir, 're', 'vw_av', re_vw)
        call Read_array(turb_data_dir, 're', "y", y)
    end subroutine

    
    function Interp(y, y_arr, val_arr)
        real :: y
        real :: y_arr(:), val_arr(:)
        real :: Interp
        integer :: i
        integer :: num
    
        num = size(y_arr)
        do i = 1, num-1
            if ((y >= y_arr(i)).AND.(y < y_arr(i + 1))) then
                Interp = val_arr(i) + (y - y_arr(i)) * (val_arr(i + 1) - val_arr(i)) / (y_arr(i + 1) - y_arr(i))
            end if
        end do
        if (y < y_arr(1)) then
            Interp = val_arr(1) - (val_arr(2) - val_arr(1)) / (y_arr(2) - y_arr(1)) * (y_arr(1) - y)
        end if
        if (y >= y_arr(num)) then
            Interp = val_arr(num) + (val_arr(num) - val_arr(num-1)) / (y_arr(num) - y_arr(num-1)) * (y - y_arr(num))
        end if
    end function
    
end module Test


    program STG_fort
    use Test    
    implicit none

    ! Variables
    integer(8) :: num_modes = 800
    real :: re_uu = 10
    real :: re_vv = 1
    real :: re_ww = 1
    real :: re_uv = 0
    real :: re_uw = 0
    real :: re_vw = 0
    real :: ls_i = 1.

    ! character(len=100) :: turb_data_dir = 'E:\Documents\tasks\others\test_sec_data\turb_data'
    ! integer(8) :: N = 69
    ! real, allocatable :: y_arr(:), re_uu_arr(:), re_vv_arr(:), re_ww_arr(:), re_uv_arr(:), re_uw_arr(:), re_vw_arr(:) 
    ! real :: y = 0.001
    
    ! ! Body of STG_fort
    ! allocate(y_arr(N))
    ! allocate(re_uu_arr(N))
    ! allocate(re_vv_arr(N))
    ! allocate(re_ww_arr(N))
    ! allocate(re_uv_arr(N))
    ! allocate(re_uw_arr(N))
    ! allocate(re_vw_arr(N))
    
    ! call Read_re_stresses(turb_data_dir, y_arr, re_uu_arr, re_vv_arr, re_ww_arr, re_uv_arr, re_uw_arr, re_vw_arr)
    ! print *, Interp(y, y_arr, re_uu_arr)
    ! print *, ''
    ! print *, y_arr
    ! print *, ''
    ! print *, re_uu_arr

    call test_Smirnov(re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, num_modes, ls_i)
    
    call test_Davidson(re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, num_modes, ls_i)
    
    
    call test_SEM( &
            re_uu, re_vv, re_ww, &
            re_uv, re_uw, re_vw, &
            1500, & 
            5., 0., 0., &
            0.5, 0.5, 0.5, &
            0.5, 0.5, 0.5, &
            0.5, 0.5, 0.5 &
        )
    
    
    end program STG_fort

