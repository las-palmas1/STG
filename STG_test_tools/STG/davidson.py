from STG.common import *
import ctypes


class STG_DavidsonData_Stationary(ctypes.Structure):
    _fields_ = [
        ('i_cnt', STG_int),
        ('j_cnt', STG_int),
        ('k_cnt', STG_int),
        ('a', ctypes.POINTER(STG_float)),
        ('b', ctypes.POINTER(STG_float)),

        ('c1', ctypes.POINTER(STG_float)),
        ('c2', ctypes.POINTER(STG_float)),
        ('c3', ctypes.POINTER(STG_float)),
        ('a11', ctypes.POINTER(STG_float)),
        ('a12', ctypes.POINTER(STG_float)),
        ('a13', ctypes.POINTER(STG_float)),
        ('a21', ctypes.POINTER(STG_float)),
        ('a22', ctypes.POINTER(STG_float)),
        ('a23', ctypes.POINTER(STG_float)),
        ('a31', ctypes.POINTER(STG_float)),
        ('a32', ctypes.POINTER(STG_float)),
        ('a33', ctypes.POINTER(STG_float)),

        ('num_modes', STG_int),
        ('energy', ctypes.POINTER(STG_float)),
        ('k_arr', ctypes.POINTER(STG_float)),
        ('u_abs', ctypes.POINTER(STG_float)),
    ]


class STG_DavidsonData_Transient(ctypes.Structure):
    _fields_ = [
        ('num_ts', STG_int),
        ('num_modes', STG_int),
        ('phi', ctypes.POINTER(STG_float)),
        ('psi', ctypes.POINTER(STG_float)),
        ('alpha', ctypes.POINTER(STG_float)),
        ('theta', ctypes.POINTER(STG_float)),

        ('u_p_prev', ctypes.POINTER(STG_float)),
        ('v_p_prev', ctypes.POINTER(STG_float)),
        ('w_p_prev', ctypes.POINTER(STG_float)),
    ]


def compute_davidson_stat_data(
        init_data: STG_InitData, num_modes: int, dissip_rate: float, visc: float, ts: float,
        stat_data: STG_DavidsonData_Stationary
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Davidson_stat_data
    func_c.argtypes = (STG_InitData, STG_int, STG_float, STG_float, STG_float,
                       ctypes.POINTER(STG_DavidsonData_Stationary))
    func_c(init_data, num_modes, dissip_rate, visc, ts, ctypes.byref(stat_data))


def free_davidson_stat_data(stat_data: STG_DavidsonData_Stationary):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_Davidson_stat_data
    func_c.argtypes = ctypes.POINTER(STG_DavidsonData_Stationary),
    func_c(ctypes.byref(stat_data))


def alloc_davidson_trans_data(
        init_data: STG_InitData, num_modes: int, num_ts_tot: int, trans_data: STG_DavidsonData_Transient
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_alloc_Davidson_trans_data
    func_c.argtypes = STG_InitData, STG_int, STG_int, ctypes.POINTER(STG_DavidsonData_Transient)
    func_c(init_data, num_modes, num_ts_tot, ctypes.byref(trans_data))


def compute_davidson_trans_data(
        stat_data: STG_DavidsonData_Stationary, num_ts_tot: int,
        trans_data: STG_DavidsonData_Transient
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Davidson_trans_data
    func_c.argtypes = STG_DavidsonData_Stationary, STG_int, ctypes.POINTER(STG_DavidsonData_Transient)
    func_c(stat_data, num_ts_tot, ctypes.byref(trans_data))


def free_davidson_trans_data(trans_data: STG_DavidsonData_Transient):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_Davidson_trans_data
    func_c.argtypes = ctypes.POINTER(STG_DavidsonData_Transient),
    func_c(ctypes.byref(trans_data))


def compute_davidson_moment_field(
        init_data: STG_InitData, stat_data: STG_DavidsonData_Stationary, trans_data: STG_DavidsonData_Transient,
        ts: float, num_ts: int, mom_field: STG_VelMomField
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Davidson_moment_field
    func_c.argtypes = (STG_InitData, STG_DavidsonData_Stationary, ctypes.POINTER(STG_DavidsonData_Transient),
                       STG_float, STG_int, ctypes.POINTER(STG_VelMomField))
    func_c(init_data, stat_data, ctypes.byref(trans_data), ts, num_ts, ctypes.byref(mom_field))


def compute_davidson_node_hist(
        init_data: STG_InitData, stat_data: STG_DavidsonData_Stationary, ts: float, num_ts_tot: int,
        trans_data: STG_DavidsonData_Transient, node_hist: STG_VelNodeHist, i: int, j: int, k: int
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Davidson_node_hist
    func_c.argtypes = (STG_InitData, STG_DavidsonData_Stationary, STG_float, STG_int,
                       ctypes.POINTER(STG_DavidsonData_Transient), ctypes.POINTER(STG_VelNodeHist),
                       STG_int, STG_int, STG_int)
    func_c(init_data, stat_data, ts, num_ts_tot, ctypes.byref(trans_data), ctypes.byref(node_hist), i, j, k)


def get_davidson_spectrum(
        delta_min, num_modes, re_uu, re_vv, re_ww, ls_i, dissip_rate, visc
):
    alpha = 1.483
    k_t = 0.5 * (re_uu + re_vv + re_ww)
    k_e = alpha * 9 * np.pi / (55 * ls_i)
    k_min = k_e / 2
    k_max = 2 * np.pi / (2 * delta_min)

    k_arr = np.linspace(k_min, k_max, num_modes)
    k_eta = dissip_rate**0.25 / visc**0.75
    u_rms = (2 / 3 * k_t)**0.5
    energy = alpha * u_rms**2 / k_e * (k_arr / k_e)**4 * \
             np.exp(-2 * (k_arr / k_eta)**2) / (1 + (k_arr / k_e)**2)**(17/6)
    return k_arr, energy
