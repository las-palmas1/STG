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

