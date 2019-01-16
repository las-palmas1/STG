import config
import ctypes
import platform
import os
import numpy as np


def search_sgt_lib(SGT_lib_name):
    if platform.system() == 'Windows':
        ext = '.dll'
    else:
        ext = '.so'
    res = ''
    for root, dirs, files in os.walk(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'projects')):
        for file in files:
            if os.path.splitext(file)[0] == SGT_lib_name and os.path.splitext(file)[1] == ext:
                res = os.path.join(root, file)
    return res


def init_rand():
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    init_rand_ctypes = ctypes.CDLL(stg_lib_fname).init_rand
    init_rand_ctypes()


def get_uniform(x1, x2, num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    get_uniform_ctypes = ctypes.CDLL(stg_lib_fname).get_uniform
    get_uniform_ctypes.restype = ctypes.POINTER(ctypes.c_float)
    get_uniform_ctypes.argtypes = ctypes.c_float, ctypes.c_float, ctypes.c_uint32

    res_p = get_uniform_ctypes(x1, x2, num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def get_normal(mu, sigma, num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    get_normal_ctypes = ctypes.CDLL(stg_lib_fname).get_normal
    get_normal_ctypes.restype = ctypes.POINTER(ctypes.c_float)
    get_normal_ctypes.argtypes = ctypes.c_float, ctypes.c_float, ctypes.c_uint32

    res_p = get_normal_ctypes(mu, sigma, num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def get_trigon(num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    get_trigon_ctypes = ctypes.CDLL(stg_lib_fname).get_trigon
    get_trigon_ctypes.restype = ctypes.POINTER(ctypes.c_float)
    get_trigon_ctypes.argtypes = ctypes.c_uint32,

    res_p = get_trigon_ctypes(num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def func(arr: list):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    test_func_ctypes = ctypes.CDLL(stg_lib_fname).test_func
    test_func_ctypes.argtypes = ctypes.POINTER(ctypes.c_float), ctypes.c_uint
    carr = (ctypes.c_float * len(arr))(*arr)
    test_func_ctypes(carr, len(arr))


class Mesh(ctypes.Structure):
    _fields_ = [
        ('x', ctypes.POINTER(ctypes.c_float)),
        ('y', ctypes.POINTER(ctypes.c_float)),
        ('z', ctypes.POINTER(ctypes.c_float))
    ]


class Scales(ctypes.Structure):
    _fields_ = [
        ('length_scale', ctypes.POINTER(ctypes.c_float)),
        ('time_scale', ctypes.POINTER(ctypes.c_float)),
    ]


class ReStress(ctypes.Structure):
    _fields_ = [
        ('re_uu', ctypes.POINTER(ctypes.c_float)),
        ('re_vv', ctypes.POINTER(ctypes.c_float)),
        ('re_ww', ctypes.POINTER(ctypes.c_float)),
        ('re_uv', ctypes.POINTER(ctypes.c_float)),
        ('re_uw', ctypes.POINTER(ctypes.c_float)),
        ('re_vw', ctypes.POINTER(ctypes.c_float)),
    ]


class InitData(ctypes.Structure):
    _fields_ = [
        ('i_cnt', ctypes.c_uint32),
        ('j_cnt', ctypes.c_uint32),
        ('k_cnt', ctypes.c_uint32),
        ('mesh', Mesh),
        ('re', ReStress),
        ('scales', Scales)
    ]


class OutData(ctypes.Structure):
    _fields_ = [
        ('time', ctypes.POINTER(ctypes.c_float)),
        ('num_ts', ctypes.c_uint32),
        ('i_cnt', ctypes.c_uint32),
        ('j_cnt', ctypes.c_uint32),
        ('k_cnt', ctypes.c_uint32),
        ('u_p', ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
        ('v_p', ctypes.POINTER(ctypes.POINTER(ctypes.c_float))),
        ('w_p', ctypes.POINTER(ctypes.POINTER(ctypes.c_float)))
    ]


class SmirnovData(ctypes.Structure):
    _fields_ = [
        ('ts', ctypes.c_float),
        ('num_ts', ctypes.c_uint32),
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


if __name__ == '__main__':
    func([1.2, 2.3, 2.5, 0.4, 2.7])
