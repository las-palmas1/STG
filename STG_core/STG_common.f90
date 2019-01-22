MODULE STG_COMMON
    
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    
    type, bind(C) :: STG_Mesh
        type(c_ptr) :: x, y, z    
    end type STG_Mesh
    
    type, bind(C) :: STG_ReStress
        type(c_ptr) :: re_uu
        type(c_ptr) :: re_vv
        type(c_ptr) :: re_ww
        type(c_ptr) :: re_uv
        type(c_ptr) :: re_uw
        type(c_ptr) :: re_vw
    end type STG_ReStress
    
    type, bind(C) :: STG_Scales
        type(c_ptr) :: length_scale
        type(c_ptr) :: time_scale
    end type STG_Scales
    
    type, bind(C) :: STG_InitData
        integer(c_long) :: i_cnt
        integer(c_long) :: j_cnt
        integer(c_long) :: k_cnt
        type(STG_Mesh) :: mesh
        type(STG_Scales) :: scales
        type(STG_ReStress) :: re
    end type STG_InitData
    
    type, bind(C) :: STG_OutDataTS
        real(c_float) :: time
        integer(c_long) :: i_cnt, j_cnt, k_cnt
        type(c_ptr) :: u_p, v_p, w_p
    end type STG_OutDataTS
    
END MODULE STG_COMMON