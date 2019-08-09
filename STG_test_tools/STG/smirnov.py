from STG.common import *
import ctypes

# TODO: необходимо разделить вспомогательные данные на стационарные и нестационарные. При этом нестационарные
#  данные должны хранится только для одного шага по времени.


class SmirnovDataTimeIndep(ctypes.Structure):
    _fields_ = [
        ('i_cnt', ctypes.c_uint32),
        ('j_cnt', ctypes.c_uint32),
        ('k_cnt', ctypes.c_uint32),
        ('c1', ctypes.POINTER(ctypes.c_float)),
        ('c2', ctypes.POINTER(ctypes.c_float)),
        ('c3', ctypes.POINTER(ctypes.c_float)),
        ('a11', ctypes.POINTER(ctypes.c_float)),
        ('a12', ctypes.POINTER(ctypes.c_float)),
        ('a13', ctypes.POINTER(ctypes.c_float)),
        ('a21', ctypes.POINTER(ctypes.c_float)),
        ('a22', ctypes.POINTER(ctypes.c_float)),
        ('a23', ctypes.POINTER(ctypes.c_float)),
        ('a31', ctypes.POINTER(ctypes.c_float)),
        ('a32', ctypes.POINTER(ctypes.c_float)),
        ('a33', ctypes.POINTER(ctypes.c_float)),

        ('num_modes', ctypes.c_uint32),
        ('k1', ctypes.POINTER(ctypes.c_float)),
        ('k2', ctypes.POINTER(ctypes.c_float)),
        ('k3', ctypes.POINTER(ctypes.c_float)),
        ('zeta1', ctypes.POINTER(ctypes.c_float)),
        ('zeta2', ctypes.POINTER(ctypes.c_float)),
        ('zeta3', ctypes.POINTER(ctypes.c_float)),
        ('xi1', ctypes.POINTER(ctypes.c_float)),
        ('xi2', ctypes.POINTER(ctypes.c_float)),
        ('xi3', ctypes.POINTER(ctypes.c_float)),
        ('omega', ctypes.POINTER(ctypes.c_float)),
        ('p1', ctypes.POINTER(ctypes.c_float)),
        ('p2', ctypes.POINTER(ctypes.c_float)),
        ('p3', ctypes.POINTER(ctypes.c_float)),
        ('q1', ctypes.POINTER(ctypes.c_float)),
        ('q2', ctypes.POINTER(ctypes.c_float)),
        ('q3', ctypes.POINTER(ctypes.c_float)),
    ]


class SmirnovDataTimeDep(ctypes.Structure):
    _fields_ = [
        ('ts', ctypes.c_float),
        ('num_ts', ctypes.c_uint32),
    ]


def compute_Smirnov_data_time_indep(init_data: InitData, num_modes: int, data: SmirnovDataTimeIndep):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).compute_Smirnov_data_time_indep
    func_c.argtypes = InitData, ctypes.c_uint32, ctypes.POINTER(SmirnovDataTimeIndep)
    func_c(init_data, num_modes, ctypes.byref(data))


def compute_Smirnov_data_time_dep(time_arr: np.ndarray, data: SmirnovDataTimeDep):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).compute_Smirnov_data_time_dep
    func_c.argtypes = ctypes.c_float, ctypes.c_uint32, ctypes.POINTER(SmirnovDataTimeDep)
    num_ts = len(time_arr) - 1
    if num_ts != 0:
        ts = time_arr[1] - time_arr[0]
    else:
        ts = 0.
    func_c(ts, num_ts, ctypes.byref(data))


def free_Smirnov_data_time_indep(data: SmirnovDataTimeIndep):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).free_Smirnov_data_time_indep
    func_c.argtypes = ctypes.POINTER(SmirnovDataTimeIndep),
    func_c(ctypes.byref(data))


def compute_Smirnov_field_ts(init_data: InitData, data_tind: SmirnovDataTimeIndep, data_tdep: SmirnovDataTimeDep,
                             out_data: OutDataTS, time_level: float):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).compute_Smirnov_field_ts
    func_c.argtypes = InitData, SmirnovDataTimeIndep, SmirnovDataTimeDep, ctypes.POINTER(OutDataTS), ctypes.c_float,
    func_c(init_data, data_tind, data_tdep, ctypes.byref(out_data), time_level)


def compute_Smirnov_field_node(init_data: InitData, data_tind: SmirnovDataTimeIndep, data_tdep: SmirnovDataTimeDep,
                               out_data: OutDataNode, i: int, j: int, k: int):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).compute_Smirnov_field_node
    func_c.argtypes = InitData, SmirnovDataTimeIndep, SmirnovDataTimeDep, ctypes.POINTER(OutDataNode), \
                      ctypes.c_uint32, ctypes.c_uint32, ctypes.c_uint32
    func_c(init_data, data_tind, data_tdep, ctypes.byref(out_data), i, j, k)


def compute_Smirnov_field(init_data: InitData, data_tind: SmirnovDataTimeIndep, data_tdep: SmirnovDataTimeDep,
                          out_data: OutData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).compute_Smirnov_field
    func_c.argtypes = InitData, SmirnovDataTimeIndep, SmirnovDataTimeDep, ctypes.POINTER(OutData)
    func_c(init_data, data_tind, data_tdep, ctypes.byref(out_data))

