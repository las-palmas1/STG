from generators.abstract import BCType, Generator, Block
from typing import List, Tuple
import numpy as np
from scipy.interpolate import interp1d

from STG.common import STG_VelNodeHist, STG_VelMomField, free_mom_field, free_node_hist, free_init_data, \
    extract_pulsations_from_mom_field, extract_pulsations_from_node_hist

from STG.smirnov import compute_smirnov_data, compute_smirnov_node_hist, compute_smirnov_moment_field, \
    free_smirnov_data, STG_SmirnovData

from generators.tools import davidson_compute_velocity_field, davidson_compute_velocity_pulsation, \
    original_sem_compute_velocity_field, original_sem_compute_pulsation

from STG.davidson import STG_DavidsonData_Transient, STG_DavidsonData_Stationary, free_davidson_trans_data, \
    free_davidson_stat_data, alloc_davidson_trans_data, compute_davidson_node_hist, compute_davidson_stat_data, \
    compute_davidson_moment_field, compute_davidson_trans_data

from STG.sem import STG_SEMData_Transient, STG_SEMData_Stationary, free_sem_stat_data, free_sem_trans_data, \
    compute_sem_stat_data, compute_sem_trans_data, compute_sem_moment_field, compute_sem_node_hist


class Lund(Generator):
    def __init__(
            self, block: Block, u_av: Tuple[float, float, float],
            re_uu: float, re_vv: float, re_ww: float,
            re_uv: float, re_uw: float, re_vw: float,
            time_arr: np.ndarray
    ):
        Generator.__init__(self, block, u_av, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, 0., 0., time_arr)

    def _compute_cholesky(self):
        """Вычисление разложения тензора рейнольдсовых напряжений по Холецкому в каждой точке."""
        self.a11 = np.sqrt(self.re_uu)
        self.a12 = 0
        self.a13 = 0
        self.a21 = self.re_uv / self.a11
        self.a22 = np.sqrt(self.re_vv - self.a21 ** 2)
        self.a23 = 0
        self.a31 = self.re_uw / self.a11
        self.a32 = (self.re_vw - self.a21 * self.a31) / self.a22
        self.a33 = np.sqrt(self.re_ww - self.a31 ** 2 - self.a32 ** 2)

    def compute_velocity_field(self, num_ts):
        u_prime = np.random.normal(0, 1, self.block.shape)
        v_prime = np.random.normal(0, 1, self.block.shape)
        w_prime = np.random.normal(0, 1, self.block.shape)
        u = self.a11 * u_prime + self.a12 * v_prime + self.a13 * w_prime
        v = self.a21 * u_prime + self.a22 * v_prime + self.a23 * w_prime
        w = self.a31 * u_prime + self.a32 * v_prime + self.a33 * w_prime
        self._vel_field = (u, v, w)

    def compute_pulsation_at_node(self):
        i = self._i_puls
        j = self._j_puls
        k = self._k_puls
        u_prime = np.random.normal(0, 1, self.time_arr.shape)
        v_prime = np.random.normal(0, 1, self.time_arr.shape)
        w_prime = np.random.normal(0, 1, self.time_arr.shape)
        u = self.a11 * u_prime + self.a12 * v_prime + self.a13 * w_prime
        v = self.a21 * u_prime + self.a22 * v_prime + self.a23 * w_prime
        w = self.a31 * u_prime + self.a32 * v_prime + self.a33 * w_prime
        self._vel_puls = (u, v, w)

    def _compute_aux_data_stationary(self):
        self._compute_cholesky()

    def _compute_aux_data_transient(self):
        pass

    def _alloc_aux_data_transient(self):
        pass

    def free_data(self):
        pass


class Smirnov(Generator):
    def __init__(self, block: Block, u_av: Tuple[float, float, float], ls_i: float, ts_i: float,
                 re_uu: float, re_vv: float, re_ww: float,
                 re_uv: float, re_uw: float, re_vw: float,
                 time_arr: np.ndarray,
                 mode_num: int = 100, ):
        self.ts_i = ts_i
        self.ls_i = ls_i
        self.mode_num = mode_num
        self._c_data = STG_SmirnovData()
        Generator.__init__(
            self, block, u_av, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
            ls_i=ls_i, ls_ux=0, ls_uy=0, ls_uz=0, ls_vx=0, ls_vy=0, ls_vz=0,
            ls_wx=0, ls_wy=0, ls_wz=0,
            ts_i=ts_i, ts_u=0, ts_v=0, ts_w=0,
            time_arr=time_arr
        )

    def _compute_aux_data_stationary(self):
        compute_smirnov_data(self._c_init_data, self.mode_num, self._c_data)

    def _compute_aux_data_transient(self):
        pass

    def compute_pulsation_at_node(self):
        i, j, k = self.get_puls_node()
        compute_smirnov_node_hist(
            init_data=self._c_init_data, data=self._c_data, ts=self._ts, num_ts=self._num_ts_tot,
            node_hist=self._c_node_hist, i=i, j=j, k=k)
        self._vel_puls = extract_pulsations_from_node_hist(self._c_node_hist)
        free_node_hist(self._c_node_hist)

    def compute_velocity_field(self, num_ts):
        compute_smirnov_moment_field(
            init_data=self._c_init_data, data=self._c_data, time=self.time_arr[num_ts], mom_field=self._c_mom_field
        )
        self._vel_field = extract_pulsations_from_mom_field(self._c_mom_field)
        free_mom_field(self._c_mom_field)

    def _get_energy_desired(self, k):
        return 16 * (2 / np.pi) ** 0.5 * k**4 * np.exp(-2 * k**2)

    def _alloc_aux_data_transient(self):
        pass

    def free_data(self):
        free_smirnov_data(self._c_data)


class Davidson(Generator):
    def __init__(
            self, block: Block, u_av: Tuple[float, float, float],
            ts_i: float, ls_i: float, dissip_rate: float, visc: float, num_modes: int,
            re_uu: float, re_vv: float, re_ww: float,
            re_uv: float, re_uw: float, re_vw: float,
            time_arr: np.ndarray,
    ):
        self.ts_i = ts_i
        self.ls_i = ls_i
        self.dissip_rate = dissip_rate
        self.visc = visc
        self.num_modes = num_modes
        self._c_stat_data = STG_DavidsonData_Stationary()
        self._c_trans_data = STG_DavidsonData_Transient()
        Generator.__init__(
            self, block, u_av, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, ls_i=ls_i,
            ls_ux=0, ls_uy=0, ls_uz=0,
            ls_vx=0, ls_vy=0, ls_vz=0,
            ls_wx=0, ls_wy=0, ls_wz=0,
            ts_i=ts_i, ts_u=0, ts_v=0, ts_w=0, time_arr=time_arr
        )

    def _alloc_aux_data_transient(self):
        alloc_davidson_trans_data(self._c_init_data, self.num_modes, self._num_ts_tot, self._c_trans_data)

    def _compute_aux_data_transient(self):
        compute_davidson_trans_data(self._c_stat_data, self._num_ts_tot, self._c_trans_data)

    def _compute_aux_data_stationary(self):
        compute_davidson_stat_data(self._c_init_data, self.num_modes, self.dissip_rate, self.visc, self._ts,
                                   self._c_stat_data)

    def compute_velocity_field(self, num_ts):
        compute_davidson_moment_field(
            self._c_init_data, self._c_stat_data, self._c_trans_data, self._ts, num_ts, self._c_mom_field
        )
        self._vel_field = extract_pulsations_from_mom_field(self._c_mom_field)
        free_mom_field(self._c_mom_field)

    def compute_pulsation_at_node(self):
        i, j, k = self.get_puls_node()
        compute_davidson_node_hist(
            self._c_init_data, self._c_stat_data, self._ts, self._num_ts_tot, self._c_trans_data, self._c_node_hist,
            i, j, k
        )
        self._vel_puls = extract_pulsations_from_node_hist(self._c_node_hist)
        free_node_hist(self._c_node_hist)

    def free_data(self):
        free_davidson_stat_data(self._c_stat_data)
        free_davidson_trans_data(self._c_trans_data)

    def _get_energy_desired(self, k):
        k_arr = np.zeros(self.num_modes)
        energy = np.zeros(self.num_modes)
        for i in range(self.num_modes):
            k_arr[i] = self._c_stat_data.k_arr[i]
            energy[i] = self._c_stat_data.energy[i]
        return float(interp1d(k_arr, energy, fill_value=0, bounds_error=False)(k))


class OriginalSEM(Generator):
    def __init__(
            self, block: Block, u_e: Tuple[float, float, float], eddies_num: int,
            ls_ux: float, ls_uy: float, ls_uz: float,
            ls_vx: float, ls_vy: float, ls_vz: float,
            ls_wx: float, ls_wy: float, ls_wz: float,
            re_uu: float, re_vv: float, re_ww: float,
            re_uv: float, re_uw: float, re_vw: float,
            time_arr: np.ndarray,
    ):
        self.eddies_num = eddies_num
        self.u_e = u_e
        self._c_stat_data = STG_SEMData_Stationary()
        self._c_trans_data = STG_SEMData_Transient()
        ls_i = (ls_ux + ls_uy + ls_uz + ls_vx + ls_vy + ls_vz + ls_wx + ls_wy + ls_wz)
        Generator.__init__(
            self, block, u_e, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
            ls_i=ls_i, ls_ux=ls_ux, ls_uy=ls_uy, ls_uz=ls_uz,
            ls_vx=ls_vx, ls_vy=ls_vy, ls_vz=ls_vz,
            ls_wx=ls_wx, ls_wy=ls_wy, ls_wz=ls_wz,
            ts_i=0., ts_u=0., ts_v=0., ts_w=0.,
            time_arr=time_arr
        )

    def _alloc_aux_data_transient(self):
        pass

    def _compute_aux_data_transient(self):
        compute_sem_trans_data(self._c_stat_data, self._ts, self._num_ts_tot, self._c_trans_data)

    def _compute_aux_data_stationary(self):
        compute_sem_stat_data(self._c_init_data, self.eddies_num, self.u_e[0], self.u_e[1],
                              self.u_e[2], self._c_stat_data)

    def compute_velocity_field(self, num_ts):
        compute_sem_moment_field(self._c_init_data, self._c_stat_data, self._c_trans_data, self._ts,
                                 self._num_ts_tot, self._c_mom_field)
        self._vel_field = extract_pulsations_from_mom_field(self._c_mom_field)
        free_mom_field(self._c_mom_field)

    def compute_pulsation_at_node(self):
        i, j, k = self.get_puls_node()
        compute_sem_node_hist(self._c_init_data, self._c_stat_data, self._c_trans_data, self._ts,
                              self._num_ts_tot, self._c_node_hist, i, j, k)
        self._vel_puls = extract_pulsations_from_node_hist(self._c_node_hist)
        free_node_hist(self._c_node_hist)

    def free_data(self):
        free_sem_stat_data(self._c_stat_data)
        free_sem_trans_data(self._c_trans_data)








