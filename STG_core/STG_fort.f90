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
        
        real, allocatable, target :: rand_data(:, :)
        real, allocatable, target :: mat_data(:)
        real, pointer :: k1(:), k2(:), k3(:), zeta1(:), zeta2(:), zeta3(:), xi1(:), xi2(:), xi3(:)
        real, pointer :: p1(:), p2(:), p3(:), q1(:), q2(:), q3(:), omega(:)
        real, pointer :: c1, c2, c3, a11, a12, a13
        real, pointer :: a21, a22, a23, a31, a32, a33
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
        ! NOTE: ссылки нужно делать на столбцы
        k1 => rand_data(:, 1)
        k2 => rand_data(:, 2)
        k3 => rand_data(:, 3)
        zeta1 => rand_data(:, 4)
        zeta2 => rand_data(:, 5)
        zeta3 => rand_data(:, 6)
        xi1 => rand_data(:, 7)
        xi2 => rand_data(:, 8)
        xi3 => rand_data(:, 9)
        p1 => rand_data(:, 10)
        p2 => rand_data(:, 11)
        p3 => rand_data(:, 12)
        q1 => rand_data(:, 13)
        q2 => rand_data(:, 14)
        q3 => rand_data(:, 15)
        omega => rand_data(:, 16)
        
        allocate(mat_data(num_mat_vars))
        c1 => mat_data(1)
        c2 => mat_data(2)
        c3 => mat_data(3)
        a11 => mat_data(4)
        a12 => mat_data(5)
        a13 => mat_data(6)
        a21 => mat_data(7)
        a22 => mat_data(8)
        a23 => mat_data(9)
        a31 => mat_data(10)
        a32 => mat_data(11)
        a33 => mat_data(12)
        
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

