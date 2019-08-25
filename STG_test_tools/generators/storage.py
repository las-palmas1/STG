from generators.abstract import BCType, Generator, Block
from typing import List, Tuple
import numpy as np
import numpy.linalg as la
from random import choices
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from STG.common import STG_VelNodeHist, STG_VelMomField, free_mom_field, free_node_hist, free_init_data, \
    extract_pulsations_from_mom_field, extract_pulsations_from_node_hist

from STG.smirnov import compute_smirnov_data, compute_smirnov_node_hist, compute_smirnov_moment_field, \
    free_smirnov_data, STG_SmirnovData

from generators.tools import davidson_compute_velocity_field, davidson_compute_velocity_pulsation, \
    original_sem_compute_velocity_field, original_sem_compute_pulsation

from STG.davidson import STG_DavidsonData_Transient, STG_DavidsonData_Stationary, free_davidson_trans_data, \
    free_davidson_stat_data, alloc_davidson_trans_data, compute_davidson_node_hist, compute_davidson_stat_data, \
    compute_davidson_moment_field, compute_davidson_trans_data


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
            self, block: Block, u_av: Tuple[float, float, float], sigma: float, eddy_num: int,
            re_uu: float, re_vv: float, re_ww: float,
            re_uv: float, re_uw: float, re_vw: float,
            time_arr: np.ndarray,
    ):
        self.sigma = sigma
        self.eddy_num = eddy_num
        self.eddy_positions_field: np.ndarray = []
        self.eddy_positions_puls: np.ndarray = []
        Generator.__init__(self, block, u_av, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, time_arr)

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

    def _compute_init_eddies_pos(self):
        self.u0 = self.u_av[0]
        self.v0 = self.u_av[1]
        self.w0 = self.u_av[2]
        self.x_min = self.block.mesh[0].min() - self.sigma
        self.x_max = self.block.mesh[0].max() + self.sigma
        self.y_min = self.block.mesh[1].min() - self.sigma
        self.y_max = self.block.mesh[1].max() + self.sigma
        self.z_min = self.block.mesh[2].min() - self.sigma
        self.z_max = self.block.mesh[2].max() + self.sigma
        self.volume = (self.x_max - self.x_min) * (self.y_max - self.y_min) * (self.z_max - self.z_min)
        # self.eddy_num = int(self.volume / self.sigma**3)
        self.x_e_init = np.random.uniform(self.x_min, self.x_max, self.eddy_num)
        self.y_e_init = np.random.uniform(self.y_min, self.y_max, self.eddy_num)
        self.z_e_init = np.random.uniform(self.z_min, self.z_max, self.eddy_num)
        self.epsilon_init = np.random.normal(0, 1, (self.eddy_num, 3))
        # self.epsilon_init = np.random.uniform(-1, 1, (self.eddy_num, 3))

    @classmethod
    def _get_line_plane_intersection(cls, x0, y0, z0, xv, yv, zv, a, b, c, d):
        k = a * xv + b * yv + c * zv
        x_res = (x0 * (b * yv + c * zv) - xv * (b * y0 + c * z0 + d))
        y_res = (y0 * (a * xv + c * zv) - yv * (a * x0 + c * z0 + d))
        z_res = (z0 * (b * yv + a * xv) - zv * (b * y0 + a * x0 + d))
        if k == 0:
            return None
        else:
            return x_res / k, y_res / k, z_res / k

    @classmethod
    def get_scalar_prod(cls, x1, y1, z1, x2, y2, z2):
        return x1 * x2 + y1 * y2 + z1 * z2

    @classmethod
    def get_in_planes(cls, x_min, x_max, y_min, y_max, z_min, z_max,
                      x_e: np.ndarray, y_e: np.ndarray, z_e: np.ndarray, u0, v0, w0):
        # коэффициенты плоскостей, ограничивающих область
        bounds_coef = (
            (1, 0, 0, -x_min), (1, 0, 0, -x_max),
            (0, 1, 0, -y_min), (0, 1, 0, -y_max),
            (0, 0, 1, -z_min), (0, 0, 1, -z_max),
        )
        # интервалы координат, в которых расположена каждая из 6 граничных граней области
        bounds = (
            ((x_min, x_min), (y_min, y_max), (z_min, z_max)),
            ((x_max, x_max), (y_min, y_max), (z_min, z_max)),
            ((x_min, x_max), (y_min, y_min), (z_min, z_max)),
            ((x_min, x_max), (y_max, y_max), (z_min, z_max)),
            ((x_min, x_max), (y_min, y_max), (z_min, z_min)),
            ((x_min, x_max), (y_min, y_max), (z_max, z_max)),
        )
        in_planes = []
        for x0, y0, z0 in zip(x_e, y_e, z_e):
            intersecs = (
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[0]),
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[1]),
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[2]),
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[3]),
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[4]),
                cls._get_line_plane_intersection(x0, y0, z0, u0, v0, w0, *bounds_coef[5]),
            )
            for intersec, bound in zip(intersecs, bounds):
                if intersec:
                    if (
                            (intersec[0] >= x_min or intersec[0] <= x_max) and
                            (intersec[1] >= y_min or intersec[1] <= y_max) and
                            (intersec[2] >= z_min or intersec[2] <= z_max)
                    ):
                        s = cls.get_scalar_prod(
                            (intersec[0] - x0), (intersec[1] - y0), (intersec[2] - z0), u0, v0, w0
                        )
                        if s < 0:
                            in_planes.append(bound)
        return in_planes

    def get_eddies_params(self, time, num_ts: int = 100):
        """
        :param time: Временной интервал, на котором нужно вычислить позиции вихрей
        :param num_ts: Число отрезков на этом интервале.

        Возвращает набор значений координат позиций вихрей и их интенсивностей в различные моменты времени.
        """
        if num_ts != 0:
            dt = time / num_ts
        else:
            dt = 0
        in_planes = self.get_in_planes(
            self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max,
            self.x_e_init, self.y_e_init, self.z_e_init, self.u0, self.v0, self.w0,
        )
        x_e = self.x_e_init
        y_e = self.y_e_init
        z_e = self.z_e_init
        epsilon = self.epsilon_init
        res = np.zeros((num_ts + 1, 6, self.eddy_num))
        for n in range(num_ts + 1):
            if n > 0:
                x_e = x_e + dt * self.u0
                y_e = y_e + dt * self.v0
                z_e = z_e + dt * self.w0
            for i in range(self.eddy_num):
                if (
                        x_e[i] < self.x_min or x_e[i] > self.x_max or
                        y_e[i] < self.y_min or y_e[i] > self.y_max or
                        z_e[i] < self.z_min or z_e[i] > self.z_max
                ):
                    x_e[i] = np.random.uniform(in_planes[i][0][0], in_planes[i][0][1])
                    y_e[i] = np.random.uniform(in_planes[i][1][0], in_planes[i][1][1])
                    z_e[i] = np.random.uniform(in_planes[i][2][0], in_planes[i][2][1])
                    epsilon[i] = np.random.normal(0, 1, 3)
                    # epsilon[i] = np.random.uniform(-1, 1, 3)

            res[n, 0, :] = x_e
            res[n, 1, :] = y_e
            res[n, 2, :] = z_e
            res[n, 3, :] = epsilon[:, 0]
            res[n, 4, :] = epsilon[:, 1]
            res[n, 5, :] = epsilon[:, 2]
        return res

    def _compute_aux_data_stationary(self):
        self._compute_init_eddies_pos()

    def _compute_aux_data_space(self):
        self._compute_cholesky()

    def _compute_aux_data_transient(self):
        self.eddy_positions_field = self.get_eddies_params(
            self.time_arr_field.max() - self.time_arr_field.min(),
            self.time_arr_field.shape[0] - 1
        )

    def compute_aux_data_transient_node(self):
        self.eddy_positions_puls = self.get_eddies_params(
            self.time_arr.max() - self.time_arr.min(),
            self.time_arr.shape[0] - 1
        )

    @classmethod
    def form_func(cls, x):
        return (np.abs(x) < 1) * (np.sqrt(1.5) * (1 - np.abs(x)))

    def compute_velocity_field(self):
        vel = original_sem_compute_velocity_field(
            positions=self.eddy_positions_field,
            x=self.block.mesh[0], y=self.block.mesh[1], z=self.block.mesh[2],
            sigma=self.sigma, volume=self.volume, eddy_num=self.eddy_num,
            a11=self.a11, a12=self.a12, a13=self.a13,
            a21=self.a21, a22=self.a22, a23=self.a23,
            a31=self.a31, a32=self.a32, a33=self.a33,
        )
        for i in range(vel.shape[0]):
            self._vel_field.append((vel[i, 0, :, :, :], vel[i, 1, :, :, :], vel[i, 2, :, :, :]))

    def compute_pulsation_at_node(self):
        i, j, k = self.get_puls_node()
        self._vel_puls = original_sem_compute_pulsation(
            positions=self.eddy_positions_puls,
            x=self.block.mesh[0][i, j, k], y=self.block.mesh[1][i, j, k], z=self.block.mesh[2][i, j, k],
            sigma=self.sigma, volume=self.volume, eddy_num=self.eddy_num,
            a11=self.a11, a12=self.a12, a13=self.a13,
            a21=self.a21, a22=self.a22, a23=self.a23,
            a31=self.a31, a32=self.a32, a33=self.a33,
        )





