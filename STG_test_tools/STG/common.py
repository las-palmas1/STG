import config
import ctypes
import platform
import os
import numpy as np
from typing import Tuple


STG_int = ctypes.c_ulong
STG_float = ctypes.c_float


def search_sgt_lib(SGT_lib_name, conf):
    if platform.system() == 'Windows':
        ext = '.dll'
    else:
        ext = '.so'
    res = ''
    for root, dirs, files in os.walk(
            os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'projects')
    ):
        for file in files:
            if os.path.splitext(file)[0] == SGT_lib_name and os.path.splitext(file)[1] == ext and \
                    root.count(conf) != 0:
                res = os.path.join(root, file)
    return res


def init_rand():
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    init_rand_ctypes = ctypes.CDLL(stg_lib_fname).STG_init_rand
    init_rand_ctypes()


def get_uniform(x1, x2, num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    get_uniform_ctypes = ctypes.CDLL(stg_lib_fname).get_uniform
    get_uniform_ctypes.restype = ctypes.POINTER(STG_float)
    get_uniform_ctypes.argtypes = STG_float, STG_float, STG_int

    res_p = get_uniform_ctypes(x1, x2, num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def get_normal(mu, sigma, num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    get_normal_ctypes = ctypes.CDLL(stg_lib_fname).get_normal
    get_normal_ctypes.restype = ctypes.POINTER(STG_float)
    get_normal_ctypes.argtypes = STG_float, STG_float, STG_int

    res_p = get_normal_ctypes(mu, sigma, num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def get_trigon(num):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    get_trigon_ctypes = ctypes.CDLL(stg_lib_fname).get_trigon
    get_trigon_ctypes.restype = ctypes.POINTER(STG_float)
    get_trigon_ctypes.argtypes = STG_int

    res_p = get_trigon_ctypes(num)
    res = np.zeros(num)
    for i in range(num):
        res[i] = res_p[i]
    return res


def func(arr: list):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    test_func_ctypes = ctypes.CDLL(stg_lib_fname).test_func
    test_func_ctypes.argtypes = ctypes.POINTER(STG_float), ctypes.c_uint
    carr = (STG_float * len(arr))(*arr)
    test_func_ctypes(carr, len(arr))


class STG_Mesh(ctypes.Structure):
    _fields_ = [
        ('x', ctypes.POINTER(STG_float)),
        ('y', ctypes.POINTER(STG_float)),
        ('z', ctypes.POINTER(STG_float))
    ]


class STG_Scales(ctypes.Structure):
    _fields_ = [
        ('ls_i', ctypes.POINTER(STG_float)),
        ('ls_ux', ctypes.POINTER(STG_float)),
        ('ls_uy', ctypes.POINTER(STG_float)),
        ('ls_uz', ctypes.POINTER(STG_float)),
        ('ls_vx', ctypes.POINTER(STG_float)),
        ('ls_vy', ctypes.POINTER(STG_float)),
        ('ls_vz', ctypes.POINTER(STG_float)),
        ('ls_wx', ctypes.POINTER(STG_float)),
        ('ls_wy', ctypes.POINTER(STG_float)),
        ('ls_wz', ctypes.POINTER(STG_float)),
        ('ts_i', ctypes.POINTER(STG_float)),
        ('ts_u', ctypes.POINTER(STG_float)),
        ('ts_v', ctypes.POINTER(STG_float)),
        ('ts_w', ctypes.POINTER(STG_float)),
    ]


class STG_ReStress(ctypes.Structure):
    _fields_ = [
        ('re_uu', ctypes.POINTER(STG_float)),
        ('re_vv', ctypes.POINTER(STG_float)),
        ('re_ww', ctypes.POINTER(STG_float)),
        ('re_uv', ctypes.POINTER(STG_float)),
        ('re_uw', ctypes.POINTER(STG_float)),
        ('re_vw', ctypes.POINTER(STG_float)),
    ]


class STG_InitData(ctypes.Structure):
    _fields_ = [
        ('i_cnt', ctypes.c_ulong),
        ('j_cnt', ctypes.c_ulong),
        ('k_cnt', ctypes.c_ulong),
        ('mesh', STG_Mesh),
        ('re', STG_ReStress),
        ('scales', STG_Scales)
    ]


class STG_VelMomField(ctypes.Structure):
    _fields_ = [
        ('time', STG_float),
        ('i_cnt', ctypes.c_ulong),
        ('j_cnt', ctypes.c_long),
        ('k_cnt', ctypes.c_long),
        ('u_p', ctypes.POINTER(STG_float)),
        ('v_p', ctypes.POINTER(STG_float)),
        ('w_p', ctypes.POINTER(STG_float))
    ]


class STG_VelNodeHist(ctypes.Structure):
    _fields_ = [
        ('time', ctypes.POINTER(STG_float)),
        ('num_ts', ctypes.c_long),
        ('i', ctypes.c_long),
        ('j', ctypes.c_long),
        ('k', ctypes.c_long),
        ('u_p', ctypes.POINTER(STG_float)),
        ('v_p', ctypes.POINTER(STG_float)),
        ('w_p', ctypes.POINTER(STG_float))
    ]


def extract_pulsations_from_mom_field(mom_field: STG_VelMomField):
    num = mom_field.i_cnt * mom_field.j_cnt * mom_field.k_cnt
    u = np.zeros(num)
    v = np.zeros(num)
    w = np.zeros(num)
    for i in range(num):
        u[i] = mom_field.u_p[i]
        v[i] = mom_field.v_p[i]
        w[i] = mom_field.w_p[i]
    u = u.reshape(mom_field.i_cnt, mom_field.j_cnt, mom_field.k_cnt)
    v = v.reshape(mom_field.i_cnt, mom_field.j_cnt, mom_field.k_cnt)
    w = w.reshape(mom_field.i_cnt, mom_field.j_cnt, mom_field.k_cnt)
    return u, v, w


def extract_pulsations_from_node_hist(node_hist: STG_VelNodeHist):
    u = np.zeros(node_hist.num_ts + 1)
    v = np.zeros(node_hist.num_ts + 1)
    w = np.zeros(node_hist.num_ts + 1)
    for i in range(node_hist.num_ts + 1):
        u[i] = node_hist.u_p[i]
        v[i] = node_hist.v_p[i]
        w[i] = node_hist.w_p[i]
    return u, v, w


def get_init_data(
        mesh: Tuple[np.ndarray, np.ndarray, np.ndarray],
        re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
        ls_i, ls_ux, ls_uy, ls_uz, ls_vx, ls_vy, ls_vz, ls_wx, ls_wy, ls_wz,
        ts_i, ts_u, ts_v, ts_w
):
    i_cnt = mesh[0].shape[0]
    j_cnt = mesh[0].shape[1]
    k_cnt = mesh[0].shape[2]
    num = i_cnt * j_cnt * k_cnt
    x = (STG_float * num)(*mesh[0].reshape(num))
    y = (STG_float * num)(*mesh[1].reshape(num))
    z = (STG_float * num)(*mesh[2].reshape(num))
    mesh_c = STG_Mesh(x, y, z)

    re_uu_arr = np.full(num, re_uu)
    re_vv_arr = np.full(num, re_vv)
    re_ww_arr = np.full(num, re_ww)
    re_uv_arr = np.full(num, re_uv)
    re_uw_arr = np.full(num, re_uw)
    re_vw_arr = np.full(num, re_vw)

    re_uu_c = (STG_float * num)(*re_uu_arr)
    re_vv_c = (STG_float * num)(*re_vv_arr)
    re_ww_c = (STG_float * num)(*re_ww_arr)
    re_uv_c = (STG_float * num)(*re_uv_arr)
    re_uw_c = (STG_float * num)(*re_uw_arr)
    re_vw_c = (STG_float * num)(*re_vw_arr)
    re_stress_c = STG_ReStress(re_uu=re_uu_c, re_vv=re_vv_c, re_ww=re_ww_c, re_uv=re_uv_c, re_uw=re_uw_c, re_vw=re_vw_c)

    ls_i_arr = np.full(num, ls_i)
    ls_ux_arr = np.full(num, ls_ux)
    ls_uy_arr = np.full(num, ls_uy)
    ls_uz_arr = np.full(num, ls_uz)
    ls_vx_arr = np.full(num, ls_vx)
    ls_vy_arr = np.full(num, ls_vy)
    ls_vz_arr = np.full(num, ls_vz)
    ls_wx_arr = np.full(num, ls_wx)
    ls_wy_arr = np.full(num, ls_wy)
    ls_wz_arr = np.full(num, ls_wz)

    ts_i_arr = np.full(num, ts_i)
    ts_u_arr = np.full(num, ts_u)
    ts_v_arr = np.full(num, ts_v)
    ts_w_arr = np.full(num, ts_w)

    ls_i_c = (STG_float * num)(*ls_i_arr)
    ls_ux_c = (STG_float * num)(*ls_ux_arr)
    ls_uy_c = (STG_float * num)(*ls_uy_arr)
    ls_uz_c = (STG_float * num)(*ls_uz_arr)
    ls_vx_c = (STG_float * num)(*ls_vx_arr)
    ls_vy_c = (STG_float * num)(*ls_vy_arr)
    ls_vz_c = (STG_float * num)(*ls_vz_arr)
    ls_wx_c = (STG_float * num)(*ls_wx_arr)
    ls_wy_c = (STG_float * num)(*ls_wy_arr)
    ls_wz_c = (STG_float * num)(*ls_wz_arr)

    ts_i_c = (STG_float * num)(*ts_i_arr)
    ts_u_c = (STG_float * num)(*ts_u_arr)
    ts_v_c = (STG_float * num)(*ts_v_arr)
    ts_w_c = (STG_float * num)(*ts_w_arr)
    scales_c = STG_Scales(
        ls_i=ls_i_c, ls_ux=ls_ux_c, ls_uy=ls_uy_c, ls_uz=ls_uz_c,
        ls_vx=ls_vx_c, ls_vy=ls_vy_c, ls_vz=ls_vz_c,
        ls_wx=ls_wx_c, ls_wy=ls_wy_c, ls_wz=ls_wz_c,
        ts_i=ts_i_c, ts_u=ts_u_c, ts_v=ts_v_c, ts_w=ts_w_c
    )

    init_data = STG_InitData(i_cnt=i_cnt, j_cnt=j_cnt, k_cnt=k_cnt, mesh=mesh_c, re=re_stress_c, scales=scales_c)
    return init_data


def free_node_hist(data: STG_VelNodeHist):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_VelNodeHist
    func_c.argtypes = ctypes.POINTER(STG_VelNodeHist),
    func_c(ctypes.byref(data))


def free_mom_field(data: STG_VelMomField):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_VelMomField
    func_c.argtypes = ctypes.POINTER(STG_VelMomField),
    func_c(ctypes.byref(data))


def free_init_data(data: STG_InitData):
    stg_lib_fname = search_sgt_lib(config.STG_lib_name, config.conf)
    func_c = ctypes.CDLL(stg_lib_fname).STG_free_InitData
    func_c.argtypes = ctypes.POINTER(STG_InitData),
    func_c(ctypes.byref(data))


if __name__ == '__main__':
    out_data = STG_VelNodeHist(i=10)
