import config
import ctypes
import platform
import os
import numpy as np
from typing import Tuple


def search_sgt_lib(SGT_lib_name):
    if platform.system() == 'Windows':
        ext = '.dll'
    else:
        ext = '.so'
    res = ''
    for root, dirs, files in os.walk(
            os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'projects')
    ):
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


class OutDataTS(ctypes.Structure):
    _fields_ = [
        ('time', ctypes.c_float),
        ('i_cnt', ctypes.c_uint32),
        ('j_cnt', ctypes.c_uint32),
        ('k_cnt', ctypes.c_uint32),
        ('u_p', ctypes.POINTER(ctypes.c_float)),
        ('v_p', ctypes.POINTER(ctypes.c_float)),
        ('w_p', ctypes.POINTER(ctypes.c_float))
    ]


class OutDataNode(ctypes.Structure):
    _fields_ = [
        ('time', ctypes.POINTER(ctypes.c_float)),
        ('num_ts', ctypes.c_uint32),
        ('i', ctypes.c_uint32),
        ('j', ctypes.c_uint32),
        ('k', ctypes.c_uint32),
        ('u_p', ctypes.POINTER(ctypes.c_float)),
        ('v_p', ctypes.POINTER(ctypes.c_float)),
        ('w_p', ctypes.POINTER(ctypes.c_float))
    ]


def extract_pulsations_ts(out_data: OutDataTS):
    num = out_data.i_cnt * out_data.j_cnt * out_data.k_cnt
    u = np.zeros(num)
    v = np.zeros(num)
    w = np.zeros(num)
    for i in range(num):
        u[i] = out_data.u_p[i]
        v[i] = out_data.v_p[i]
        w[i] = out_data.w_p[i]
    u = u.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt)
    v = v.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt)
    w = w.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt)
    return u, v, w


def extract_pulsations_node(out_data: OutDataNode):
    u = np.zeros(out_data.num_ts)
    v = np.zeros(out_data.num_ts)
    w = np.zeros(out_data.num_ts)
    for i in range(out_data.num_ts + 1):
        u[i] = out_data.u_p[i]
        v[i] = out_data.v_p[i]
        w[i] = out_data.w_p[i]
    return u, v, w


def extract_pulsations(out_data: OutData):
    num = out_data.i_cnt * out_data.j_cnt * out_data.k_cnt
    vel = []
    u = np.zeros(num)
    v = np.zeros(num)
    w = np.zeros(num)
    for it in range(out_data.num_ts + 1):
        for i in range(num):
            u[i] = out_data.u_p[it][i]
            v[i] = out_data.v_p[it][i]
            w[i] = out_data.w_p[it][i]
        vel.append((
            u.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt),
            v.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt),
            w.reshape(out_data.i_cnt, out_data.j_cnt, out_data.k_cnt),
        ))
    return vel


def get_init_data(mesh: Tuple[np.ndarray, np.ndarray, np.ndarray],
                  re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, l_t, tau_t):
    i_cnt = mesh[0].shape[0]
    j_cnt = mesh[0].shape[1]
    k_cnt = mesh[0].shape[2]
    num = i_cnt * j_cnt * k_cnt
    x = (ctypes.c_float * num)(*mesh[0].reshape(num))
    y = (ctypes.c_float * num)(*mesh[1].reshape(num))
    z = (ctypes.c_float * num)(*mesh[2].reshape(num))
    mesh_c = Mesh(x=x, y=y, z=z)

    re_uu_arr = np.full(num, re_uu)
    re_vv_arr = np.full(num, re_vv)
    re_ww_arr = np.full(num, re_ww)
    re_uv_arr = np.full(num, re_uv)
    re_uw_arr = np.full(num, re_uw)
    re_vw_arr = np.full(num, re_vw)

    re_uu_c = (ctypes.c_float * num)(*re_uu_arr)
    re_vv_c = (ctypes.c_float * num)(*re_vv_arr)
    re_ww_c = (ctypes.c_float * num)(*re_ww_arr)
    re_uv_c = (ctypes.c_float * num)(*re_uv_arr)
    re_uw_c = (ctypes.c_float * num)(*re_uw_arr)
    re_vw_c = (ctypes.c_float * num)(*re_vw_arr)
    re_stress_c = ReStress(re_uu=re_uu_c, re_vv=re_vv_c, re_ww=re_ww_c, re_uv=re_uv_c, re_uw=re_uw_c, re_vw=re_vw_c)

    l_t_arr = np.full(num, l_t)
    tau_t_arr = np.full(num, tau_t)
    l_t_c = (ctypes.c_float * num)(*l_t_arr)
    tau_t_c = (ctypes.c_float * num)(*tau_t_arr)
    scales_c = Scales(length_scale=l_t_c, time_scale=tau_t_c)

    init_data = InitData(i_cnt=i_cnt, j_cnt=j_cnt, k_cnt=k_cnt, mesh=mesh_c, re=re_stress_c, scales=scales_c)
    return init_data


def free_out_data(data: OutData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).free_OutData
    func_c.argtypes = ctypes.POINTER(OutData),
    func_c(ctypes.byref(data))


def free_out_data_node(data: OutDataNode):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).free_OutDataNode
    func_c.argtypes = ctypes.POINTER(OutDataNode),
    func_c(ctypes.byref(data))


def free_out_data_ts(data: OutDataTS):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).free_OutDataTS
    func_c.argtypes = ctypes.POINTER(OutDataTS),
    func_c(ctypes.byref(data))


def free_init_data(data: InitData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name)
    func_c = ctypes.CDLL(stg_lib_fname).free_InitData
    func_c.argtypes = ctypes.POINTER(InitData),
    func_c(ctypes.byref(data))


if __name__ == '__main__':
    out_data = OutDataNode(i=10)
