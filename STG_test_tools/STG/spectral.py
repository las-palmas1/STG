from STG.common import *
import ctypes


class STG_SpectralData(ctypes.Structure):
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
        ('u_abs', ctypes.POINTER(STG_float)),

        ('phi', ctypes.POINTER(STG_float)),
        ('psi', ctypes.POINTER(STG_float)),
        ('alpha', ctypes.POINTER(STG_float)),
        ('theta', ctypes.POINTER(STG_float)),
        ('omega', ctypes.POINTER(STG_float))
    ]


def compute_spectral_data(init_data: STG_InitData, num_modes: int, data: STG_SpectralData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Spectral_data
    func_c.argtypes = STG_InitData, STG_int, ctypes.POINTER(STG_SpectralData)
    func_c(init_data, num_modes, ctypes.byref(data))


def free_spectral_data(data: STG_SpectralData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_Spectral_data
    func_c.argtypes = ctypes.POINTER(STG_SpectralData),
    func_c(ctypes.byref(data))


def compute_spectral_moment_field(
        init_data: STG_InitData, data: STG_SpectralData, ts: float, num_ts, mom_field: STG_VelMomField
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Spectral_moment_field
    func_c.argtypes = STG_InitData, STG_SpectralData, STG_float, STG_int, ctypes.POINTER(STG_VelMomField)
    func_c(init_data, data, ts, num_ts, ctypes.byref(mom_field))


def compute_spectral_node_hist(
        init_data: STG_InitData, data: STG_SpectralData, ts: float, num_ts_tot: int, node_hist: STG_VelNodeHist,
        i: int, j: int, k: int
):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).STG_compute_Spectral_node_hist
    func_c.argtypes = STG_InitData, STG_SpectralData, STG_float, STG_int, ctypes.POINTER(STG_VelNodeHist), \
                      STG_int, STG_int, STG_int
    func_c(init_data, data, ts, num_ts_tot, ctypes.byref(node_hist), i, j, k)
