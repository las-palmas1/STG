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

    program STG_fort
    
    !use STG_COMMON_STRUCT
    use STG_COMMON
    !use STG_SMIRNOV_STRUCT
    use STG_SMIRNOV
    
    implicit none

    ! Variables

    integer :: num_modes = 20
    integer :: num_rand_vars = 16
    integer :: num_mat_vars = 12
    real, allocatable, target :: rand_data(:, :)
    real, allocatable, target :: mat_data(:)
    real, pointer :: k1(:), k2(:), k3(:), zeta1(:), zeta2(:), zeta3(:), xi1(:), xi2(:), xi3(:)
    real, pointer :: p1(:), p2(:), p3(:), q1(:), q2(:), q3(:), omega(:)
    
    ! Body of STG_fort
    
    allocate(rand_data(num_rand_vars, num_modes))
    k1 => rand_data(1, :)
    k2 => rand_data(2, :)
    k3 => rand_data(3, :)
    zeta1 => rand_data(4, :)
    zeta2 => rand_data(5, :)
    zeta3 => rand_data(6, :)
    xi1 => rand_data(7, :)
    xi2 => rand_data(8, :)
    xi3 => rand_data(9, :)
    p1 => rand_data(10, :)
    p2 => rand_data(11, :)
    p3 => rand_data(12, :)
    q1 => rand_data(13, :)
    q2 => rand_data(14, :)
    q3 => rand_data(15, :)
    omega => rand_data(16, :)

    call STG_init_rand()
    call STG_compute_Smirnov_random_data( &
        num_modes, c_loc(omega(1)), &
        c_loc(k1(1)), c_loc(k2(1)), c_loc(k3(1)), &
        c_loc(zeta1(1)), c_loc(zeta2(1)), c_loc(zeta3(1)), &
        c_loc(xi1(1)), c_loc(xi2(1)), c_loc(xi3(1)), &
        c_loc(p1(1)), c_loc(p2(1)), c_loc(p3(1)), &
        c_loc(q1(1)), c_loc(q2(1)), c_loc(q3(1)) &
    )
    
    print *, k1(1), ' ',  k1(2)
    print *, p1(1), ' ',  p1(2)
    print *, q1(1), ' ',  q2(1)
    
    end program STG_fort

