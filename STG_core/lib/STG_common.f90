module MSTG_Common
    use, intrinsic :: ISO_C_BINDING
    implicit none
    
    interface 
        subroutine STG_init_rand() & 
        bind(C, name='STG_init_rand')
        end subroutine STG_init_rand
    end interface
    

end module MSTG_Common