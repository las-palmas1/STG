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
    
    ! TODO: Думаю, что для случая задания ГУ нужно создать свои структуры для хранения вспомогательных данных SmirnovData и начальных данных.
    ! Не нужно хранить в них координаты, а также рейнольдсовы напряжения и масштабы в виде массивов для каждого узла, 
    ! так как координаты проще извлечь непосредственно в процедуре расчета ГУ, а масштабы и рейнодсовы напряжения при расчете ГУ не зависят от координаты.
    ! Текущая функция для расчета пульсаций на шаге по времени уже понадобится
    
    use STG_COMMON_STRUCT
    use STG_COMMON_FUNC
    use STG_SMIRNOV_STRUCT
    use STG_SMIRNOV_FUNC
    
    implicit none

    ! Variables
    integer :: I
    integer, parameter :: size = 5
    integer :: num_modes = 20
    integer, parameter :: sum_size = size * size * size
    
    type(STG_InitData_C) :: init_data
    type(STG_SmirnovDataTimeIndep_C) :: data_tind
    
    real :: X(sum_size) = (/ (1, I=1, sum_size) /)
    real :: Y(sum_size) = (/ (2, I=1, sum_size) /)
    real :: Z(sum_size) = (/ (3, I=1, sum_size) /)
    real :: l_t(sum_size) = (/ (1, I=1, sum_size) /) 
    real :: tau_t(sum_size) = (/ (2, I=1, sum_size) /)
    real :: re_uu(sum_size) = (/ (3, I=1, sum_size) /)
    real :: re_vv(sum_size) = (/ (4, I=1, sum_size) /)
    real :: re_ww(sum_size) = (/ (5, I=1, sum_size) /)
    real :: re_uv(sum_size) = (/ (0, I=1, sum_size) /)
    real :: re_uw(sum_size) = (/ (0, I=1, sum_size) /)
    real :: re_vw(sum_size) = (/ (0, I=1, sum_size) /)
    
    real, pointer :: py(:)
    real, pointer :: p_re_uu(:), p_p1(:)
    type(c_ptr) :: cp
    
    ! Body of STG_fort
    
    init_data%i_cnt = size
    init_data%j_cnt = size
    init_data%k_cnt = size
    
    init_data%mesh%x = C_LOC(X(1))
    init_data%mesh%y = C_LOC(Y(1))
    init_data%mesh%z = C_LOC(Z(1))
    
    init_data%re%re_uu = C_LOC(re_uu(1))
    init_data%re%re_vv = C_LOC(re_vv(1))
    init_data%re%re_ww = C_LOC(re_ww(1))
    init_data%re%re_uv = C_LOC(re_uv(1))
    init_data%re%re_uw = C_LOC(re_uw(1))
    init_data%re%re_vw = C_LOC(re_vw(1))
    
    init_data%scales%length_scale = C_LOC(l_t(1))
    init_data%scales%time_scale = C_LOC(tau_t(1))
    
    call STG_init_rand()
    call STG_compute_Smirnov_data_time_indep_C(init_data, num_modes, C_LOC(data_tind))
    
    call c_f_pointer(init_data%mesh%y, py, shape=[sum_size])
    call c_f_pointer(init_data%re%re_uu, p_re_uu, shape=[sum_size])
    call c_f_pointer(data_tind%p1, p_p1, shape=[num_modes])
    
    print *, X(1)
    print *, py(1)
    print *, p_re_uu(1)
    print *, p_p1(1), ' ', p_p1(2)
    
    end program STG_fort

