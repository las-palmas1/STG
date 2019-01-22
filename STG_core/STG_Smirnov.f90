module STG_SMIRNOV
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    type, bind(C) :: STG_SmirnovDataTimeIndep
        integer(c_long) :: i_cnt, j_cnt, k_cnt
        type(c_ptr) :: c1, c2, c3
        type(c_ptr) :: a11, a12, a13
        type(c_ptr) :: a21, a22, a23
        type(c_ptr) :: a31, a32, a33
        
        integer(c_long) :: num_modes
        type(c_ptr) :: k1, k2, k3
        type(c_ptr) :: zeta1, zeta2, zeta3
        type(c_ptr) :: xi1, xi2, xi3
        type(c_ptr) :: omega
        type(c_ptr) :: p1, p2, p3
        type(c_ptr) :: q1, q2, q3
    end type STG_SmirnovDataTimeIndep
    
    type, bind(C) :: STG_SmirnovDataTimeDep
        real(c_float) :: ts
        integer(c_long) :: num_ts
    end type STG_SmirnovDataTimeDep
    
    
    interface 
        subroutine STG_compute_Smirnov_data_time_indep(init_data, num_modes, data_tind) & 
        bind(C, name='STG_compute_Smirnov_data_time_indep')
            use STG_COMMON
            type(STG_InitData), intent(in), value :: init_data
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: data_tind
        end subroutine STG_compute_Smirnov_data_time_indep
    end interface 
    
    
end module STG_SMIRNOV