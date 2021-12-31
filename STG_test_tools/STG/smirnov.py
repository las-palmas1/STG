from STG.common import *
import ctypes


class STG_SmirnovData(ctypes.Structure):
    _fields_ = [
        ('i_cnt', STG_int),
        ('j_cnt', STG_int),
        ('k_cnt', STG_int),
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
        ('k1', ctypes.POINTER(STG_float)),
        ('k2', ctypes.POINTER(STG_float)),
        ('k3', ctypes.POINTER(STG_float)),
        ('zeta1', ctypes.POINTER(STG_float)),
        ('zeta2', ctypes.POINTER(STG_float)),
        ('zeta3', ctypes.POINTER(STG_float)),
        ('xi1', ctypes.POINTER(STG_float)),
        ('xi2', ctypes.POINTER(STG_float)),
        ('xi3', ctypes.POINTER(STG_float)),
        ('omega', ctypes.POINTER(STG_float)),
        ('p1', ctypes.POINTER(STG_float)),
        ('p2', ctypes.POINTER(STG_float)),
        ('p3', ctypes.POINTER(STG_float)),
        ('q1', ctypes.POINTER(STG_float)),
        ('q2', ctypes.POINTER(STG_float)),
        ('q3', ctypes.POINTER(STG_float)),
    ]


def compute_smirnov_data(init_data: STG_InitData, num_modes: int, data: STG_SmirnovData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Smirnov_data
    func_c.argtypes = STG_InitData, STG_int, ctypes.POINTER(STG_SmirnovData)
    func_c(init_data, num_modes, ctypes.byref(data))


def free_smirnov_data(data: STG_SmirnovData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_Smirnov_data
    func_c.argtypes = ctypes.POINTER(STG_SmirnovData),
    func_c(ctypes.byref(data))


def compute_smirnov_moment_field(
        init_data: STG_InitData, data: STG_SmirnovData, time: float,
        mom_field: STG_VelMomField
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Smirnov_moment_field
    func_c.argtypes = STG_InitData, STG_SmirnovData, STG_float, ctypes.POINTER(STG_VelMomField)
    func_c(init_data, data, time, ctypes.byref(mom_field))


def compute_smirnov_node_hist(
        init_data: STG_InitData, data: STG_SmirnovData, ts: float, num_ts: int,
        node_hist: STG_VelNodeHist, i: int, j: int, k: int
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Smirnov_node_hist
    func_c.argtypes = STG_InitData, STG_SmirnovData, STG_float, STG_int, ctypes.POINTER(STG_VelNodeHist), \
                      STG_int, STG_int, STG_int
    func_c(init_data, data, ts, num_ts, ctypes.byref(node_hist), i, j, k)

