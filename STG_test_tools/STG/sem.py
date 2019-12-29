from STG.common import *
import ctypes


class Vector(ctypes.Structure):
    _fields_ = [
        ('x', STG_float),
        ('y', STG_float),
        ('z', STG_float),
    ]


class Limits(ctypes.Structure):
    _fields_ = [
        ('x_min', STG_float),
        ('x_max', STG_float),
        ('y_min', STG_float),
        ('y_max', STG_float),
        ('z_min', STG_float),
        ('z_max', STG_float),
    ]


class STG_SEMData_Stationary(ctypes.Structure):
    _fields_ = [
        ('i_cnt', STG_int),
        ('j_cnt', STG_int),
        ('k_cnt', STG_int),

        ('a11', ctypes.POINTER(STG_float)),
        ('a12', ctypes.POINTER(STG_float)),
        ('a13', ctypes.POINTER(STG_float)),
        ('a21', ctypes.POINTER(STG_float)),
        ('a22', ctypes.POINTER(STG_float)),
        ('a23', ctypes.POINTER(STG_float)),
        ('a31', ctypes.POINTER(STG_float)),
        ('a32', ctypes.POINTER(STG_float)),
        ('a33', ctypes.POINTER(STG_float)),

        # NOTE: порядок должен быть такой же, как в СИ
        ('num_eddies', STG_int),
        ('vol_lims', Limits),
        ('eddies_vel', Vector),
        ('eddies_pos_init', ctypes.POINTER(Vector)),
        ('eddies_int_init', ctypes.POINTER(Vector)),
        ('in_planes_lims', ctypes.POINTER(Limits)),
    ]


class STG_SEMData_Transient(ctypes.Structure):
    _fields_ = [
        ('num_ts', STG_int),
        ('ts', STG_float),
        ('num_eddies', STG_int),
        ('eddies_int', ctypes.POINTER(Vector)),
        ('eddies_pos', ctypes.POINTER(Vector)),
    ]


def compute_sem_stat_data(
        init_data: STG_InitData, num_eddies: int, u_e: float, v_e: float, w_e: float,
        stat_data: STG_SEMData_Stationary
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_SEM_stat_data
    func_c.argtypes = STG_InitData, STG_int, Vector, ctypes.POINTER(STG_SEMData_Stationary)
    vel = Vector(x=u_e, y=v_e, z=w_e)
    func_c(init_data, num_eddies, vel, ctypes.byref(stat_data))


def free_sem_stat_data(stat_data: STG_SEMData_Stationary):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_SEM_stat_data
    func_c.argtypes = ctypes.POINTER(STG_SEMData_Stationary),
    func_c(stat_data)


def compute_sem_trans_data(
        stat_data: STG_SEMData_Stationary, ts: float, num_ts: int, trans_data: STG_SEMData_Transient
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_SEM_trans_data
    func_c.argtypes = STG_SEMData_Stationary, STG_float, STG_int, ctypes.POINTER(STG_SEMData_Transient)
    func_c(stat_data, ts, num_ts, ctypes.byref(trans_data))


def free_sem_trans_data(trans_data: STG_SEMData_Transient):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_SEM_trans_data
    func_c.argtypes = ctypes.POINTER(STG_SEMData_Transient),
    func_c(ctypes.byref(trans_data))


def compute_sem_moment_field(
        init_data: STG_InitData, stat_data: STG_SEMData_Stationary, trans_data: STG_SEMData_Transient,
        ts: float, num_ts: int, mom_field: STG_VelMomField
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_SEM_moment_field
    func_c.argtypes = (STG_InitData, STG_SEMData_Stationary, STG_SEMData_Transient,
                       STG_float, STG_int, ctypes.POINTER(STG_VelMomField))
    func_c(init_data, stat_data, trans_data, ts, num_ts, ctypes.byref(mom_field))


def compute_sem_node_hist(
        init_data: STG_InitData, stat_data: STG_SEMData_Stationary,
        trans_data: STG_SEMData_Transient, ts: float, num_ts: int, node_hist: STG_VelNodeHist,
        i: int, j: int, k: int
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_SEM_node_hist
    func_c.argtypes = (STG_InitData, STG_SEMData_Stationary, STG_SEMData_Transient, STG_float,
                       STG_int, ctypes.POINTER(STG_VelNodeHist), STG_int, STG_int, STG_int)
    func_c(init_data, stat_data, trans_data, ts, num_ts, ctypes.byref(node_hist), i, j, k)
