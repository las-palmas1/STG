module STG_COMMON_STRUCT
    
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    
    type, bind(C) :: STG_Mesh_C
        type(c_ptr) :: x, y, z    
    end type STG_Mesh_C
    
    type  :: STG_Mesh
        real, allocatable :: x(:), y(:), z(:)
    end type STG_Mesh
    
    type, bind(C) :: STG_ReStress_C
        type(c_ptr) :: re_uu
        type(c_ptr) :: re_vv
        type(c_ptr) :: re_ww
        type(c_ptr) :: re_uv
        type(c_ptr) :: re_uw
        type(c_ptr) :: re_vw
    end type STG_ReStress_C
    
    type :: STG_ReStress
        real, allocatable :: re_uu(:), re_vv(:), re_ww(:)
        real, allocatable :: re_uv(:), re_uw(:), re_vw(:)
    end type STG_ReStress
    
    type, bind(C) :: STG_Scales_C
        type(c_ptr) :: length_scale
        type(c_ptr) :: time_scale
    end type STG_Scales_C
    
    type :: STG_Scales
        real, allocatable :: length_scale, time_scale
    end type STG_scales
    
    type, bind(C) :: STG_InitData_C
        integer(c_long) :: i_cnt
        integer(c_long) :: j_cnt
        integer(c_long) :: k_cnt
        type(STG_Mesh_C) :: mesh
        type(STG_Scales_C) :: scales
        type(STG_ReStress_C) :: re
    end type STG_InitData_C
    
    type :: STG_InitData
        integer(8) :: i_cnt
        integer(8) :: j_cnt
        integer(8) :: k_cnt
        type(STG_Mesh) :: mesh
        type(STG_Scales) :: scales
        type(STG_ReStress) :: re
    end type STG_InitData
    
    type, bind(C) :: STG_OutDataTS_C
        real(c_float) :: time
        integer(c_long) :: i_cnt, j_cnt, k_cnt
        type(c_ptr) :: u_p, v_p, w_p
    end type STG_OutDataTS_C
    
    type :: STG_OutDataTS
        real :: time
        integer(8) :: i_cnt, j_cnt, k_cnt
        real, allocatable :: u_p, v_p, w_p
    end type STG_OutDataTS
    
end module STG_COMMON_STRUCT
    
    
module STG_COMMON_FUNC
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    interface 
        subroutine STG_init_rand() & 
        bind(C, name='STG_init_rand')
        end subroutine STG_init_rand
    end interface
    
   ! contains 

end module STG_COMMON_FUNC