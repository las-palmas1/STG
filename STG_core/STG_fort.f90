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
    
    subroutine test_Smirnov( &
        re_uu, re_vv, re_ww, &
        re_uv, re_uw, re_vw, &
        num_modes &
    )
        use STG_COMMON
        use STG_SMIRNOV
        implicit none
        
        real :: re_uu, re_vv, re_ww
        real :: re_uv, re_uw, re_vw
        integer(8) :: num_modes
        
        integer :: num_rand_vars, num_mat_vars
        real :: time, x, y, z, length_scale, time_scale
        
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
        length_scale = 1
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
    end subroutine test_Smirnov
    
    
    program STG_fort
    
    implicit none

    ! Variables
    integer(8) :: num_modes = 50
    real :: re_uu = 1
    real :: re_vv = 1
    real :: re_ww = 1
    real :: re_uv = 2
    real :: re_uw = 1
    real :: re_vw = -3

    ! Body of STG_fort
    
    call test_Smirnov(re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, num_modes)
    
    end program STG_fort

