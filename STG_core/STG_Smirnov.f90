module STG_SMIRNOV_STRUCT
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    type, bind(C) :: STG_SmirnovDataTimeIndep_C
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
    end type STG_SmirnovDataTimeIndep_C
    
    type, bind(C) :: STG_SmirnovDataTimeDep_C
        real(c_float) :: ts
        integer(c_long) :: num_ts
    end type STG_SmirnovDataTimeDep_C
    
end module STG_SMIRNOV_STRUCT

    
module STG_SMIRNOV_FUNC
    
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    interface 
        subroutine STG_compute_Smirnov_data_time_indep_C(init_data, num_modes, p_data_tind) & 
        bind(C, name='STG_compute_Smirnov_data_time_indep')
            use STG_COMMON_STRUCT
            type(STG_InitData_C), intent(in), value :: init_data
            integer(c_long), intent(in), value :: num_modes
            type(c_ptr), value :: p_data_tind
        end subroutine STG_compute_Smirnov_data_time_indep_C
    end interface 
    
    interface 
        subroutine STG_compute_Smirnov_data_time_dep_C(ts, num_ts, p_data_tdep) & 
        bind(C, name='STG_compute_Smirnov_data_time_dep')
            use STG_COMMON_STRUCT
            real(c_float), intent(in), value :: ts
            integer(c_long), intent(in), value :: num_ts
            type(c_ptr), value :: p_data_tdep
        end subroutine STG_compute_Smirnov_data_time_dep_C
    end interface 
    
    interface 
        subroutine STG_free_Smirnov_data_time_indep_C(p_data_tind) & 
        bind(C, name='STG_free_Smirnov_data_time_indep')
            use STG_COMMON_STRUCT
            type(c_ptr), value :: p_data_tind
        end subroutine STG_free_Smirnov_data_time_indep_C
    end interface 
    
    interface 
        subroutine STG_compute_Smirnov_field_ts(init_data, data_tind, data_tdep, p_out_data, time_level) & 
        bind(C, name='STG_compute_Smirnov_field_ts')
            use STG_SMIRNOV_STRUCT
            use STG_COMMON_STRUCT
            type(STG_InitData_C), intent(in), value :: init_data
            type(STG_SmirnovDataTimeIndep_C), intent(in), value :: data_tind
            type(STG_SmirnovDataTimeDep_C), intent(in), value :: data_tdep
            type(c_ptr), value :: p_out_data
            integer(c_long), intent(in), value :: time_level
        end subroutine STG_compute_Smirnov_field_ts
    end interface 
    
end module STG_SMIRNOV_FUNC