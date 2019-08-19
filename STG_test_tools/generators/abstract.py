import numpy as np
from abc import ABCMeta, abstractmethod
from typing import List, Tuple
import enum
from STG.common import get_init_data


class BCType(enum.Enum):
    NotWall = 0
    Wall = 1


class Block:
    """
    Пока предполагается, что в данном классе будет содержаться информация об ортогональном трехмерном блоке,
    включая данные о типе его границ (стенка или не стенка). Сетка может быть неравномерной по всем направлениям.
    Также должны расчитываться шаги сетки по всем направлениям в каждом узле и расстояние до стенки в каждом узле.
    """
    def __init__(self, shape: tuple, mesh: Tuple[np.ndarray, np.ndarray, np.ndarray],
                 bc: List[Tuple[BCType, BCType]]):
        self.shape = shape
        self.mesh = mesh
        self.bc = bc
        self.dim = len(shape)
        self._check_dim()

    def _check_dim(self):
        if len(self.bc) != self.dim:
            raise Exception("Incorrect parameters of block")
        for i in self.mesh:
            if i.shape != self.shape:
                raise Exception("Incorrect parameters of block")


class Generator(metaclass=ABCMeta):
    """
    Базовый класс, в котором должен быть определен метод для расчета поля однородной анизотропной турбулентности
    в заданные моменты времени на 2-D или 3-D сетке. Однородность означает, что параметры турбулентности будут
    постоянны во всех точках области (масштабы, тензор напряжений Рейнольдса и др.),
    а анизотропия будет означать, что матрица тензора напряжений Рейнольдса может быть не только
    пропорцианальной единичной.
    """
    def __init__(
            self, block: Block, u_av: Tuple[float, float, float],
            re_uu: float, re_vv: float, re_ww: float,
            re_uv: float, re_uw: float, re_vw: float,
            ls_i: float, ls_ux: float, ls_uy: float, ls_uz: float,
            ls_vx: float, ls_vy: float, ls_vz: float,
            ls_wx: float, ls_wy: float, ls_wz: float,
            ts_i: float, ts_u: float, ts_v: float, ts_w: float,
            time_arr: np.ndarray,
            ):
        """
        :param block: Блок сетки, на которой нужно генерировать пульсации.
        :param u_av: Кортеж из трех значений составляющих осредненной скорости.
        :param re_uu: Осредненное произведение vx*vx.
        :param re_vv: Осредненное произведение vy*vy.
        :param re_ww: Осредненное произведение vz*vz.
        :param re_uv: Осредненное произведение vx*vy.
        :param re_uw: Осредненное произведение vx*vz.
        :param re_vw: Осредненное произведение vy*vz.
        :param ls_i: Интегральный линейный масштаб турбулентности.
        :param ts_i: Интегральный временной масштаб турбулентности.
        :param time_arr: Моменты времени, для которых производится вычисление поля скоростей.
        """
        self.block = block
        self.u_av = u_av
        self.re_uu = re_uu
        self.re_vv = re_vv
        self.re_ww = re_ww
        self.re_uv = re_uv
        self.re_uw = re_uw
        self.re_vw = re_vw
        self.ls_i = ls_i
        self.ls_ux = ls_ux
        self.ls_uy = ls_uy
        self.ls_uz = ls_uz
        self.ls_vx = ls_vx
        self.ls_vy = ls_vy
        self.ls_vz = ls_vz
        self.ls_wx = ls_wx
        self.ls_wy = ls_wy
        self.ls_wz = ls_wz
        self.ts_i = ts_i
        self._time_arr = time_arr
        self._i_puls = 0
        self._j_puls = 0
        self._k_puls = 0
        self._c_mom_filed = None
        self._c_node_hist = None
        self._vel_field: Tuple[np.ndarray, np.ndarray, np.ndarray] = ()
        self._vel_puls: Tuple[np.ndarray, np.ndarray, np.ndarray] = (
            np.zeros(time_arr.shape), np.zeros(time_arr.shape), np.zeros(time_arr.shape)
        )
        self._c_init_data = get_init_data(
            self.block.mesh, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
            ls_i=ls_i, ls_ux=ls_ux, ls_uy=ls_uy, ls_uz=ls_uz, ls_vx=ls_vx, ls_vy=ls_vy, ls_vz=ls_vz,
            ls_wx=ls_wx, ls_wy=ls_wy, ls_wz=ls_wz, ts_i=ts_i, ts_u=ts_u, ts_v=ts_v, ts_w=ts_w
        )
        self._ts = time_arr[1] - time_arr[0]
        self._num_ts = time_arr.shape[0] - 1
        self._compute_aux_data_stationary()

    @abstractmethod
    def free_data(self):
        pass

    def _get_energy_desired(self, k) -> float:
        """Для спектральных методов следует переопределить. По умолчанию возвращает ноль."""
        return 0.

    def get_desired_spectrum(self, k_arr: np.ndarray) -> np.ndarray:
        E_arr = np.array([self._get_energy_desired(k) for k in k_arr])
        return E_arr

    @classmethod
    def get_divergence(cls, vel: Tuple[np.ndarray, np.ndarray, np.ndarray],
                       mesh: Tuple[np.ndarray, np.ndarray, np.ndarray],
                       shape: tuple) -> np.ndarray:
        res = np.zeros(shape)
        for i in range(shape[0] - 1):
            for j in range(shape[1] - 1):
                if shape[2] > 1:
                    for k in range(shape[2] - 1):
                        dvx_dx = (vel[0][i + 1, j, k] - vel[0][i, j, k]) / (mesh[0][i + 1, j, k] - mesh[0][i, j, k])
                        dvy_dy = (vel[1][i, j + 1, k] - vel[1][i, j, k]) / (mesh[1][i, j + 1, k] - mesh[1][i, j, k])
                        dvz_dz = (vel[2][i, j, k + 1] - vel[2][i, j, k]) / (mesh[2][i, j, k + 1] - mesh[2][i, j, k])
                        res[i, j, k] = dvx_dx + dvy_dy + dvz_dz
                else:
                    dvx_dx = (vel[0][i + 1, j, 0] - vel[0][i, j, 0]) / (mesh[0][i + 1, j, 0] - mesh[0][i, j, 0])
                    dvy_dy = (vel[1][i, j + 1, 0] - vel[1][i, j, 0]) / (mesh[1][i, j + 1, 0] - mesh[1][i, j, 0])
                    dvz_dz = 0
                    res[i, j, 0] = dvx_dx + dvy_dy + dvz_dz
        res[shape[0]-1, 0: shape[1]-1, 0: shape[2]-1] = res[shape[0]-2, 0: shape[1]-1, 0: shape[2]-1]
        res[0: shape[0], shape[1]-1, 0: shape[2]-1] = res[0: shape[0], shape[1]-2, 0: shape[2]-1]
        if shape[2] > 1:
            res[0: shape[0], 0: shape[1], shape[2]-1] = res[0: shape[0], 0: shape[1], shape[2]-2]
        return res

    @abstractmethod
    def _compute_aux_data_stationary(self):
        """
        Вычисление вспомогательных данных, независимых от времени.
        """
        pass

    @abstractmethod
    def compute_aux_data_transient_puls(self):
        """
        Вычисление вспомогательных данных зависимых от времени
        для моментов времени, в которые вычисляются скорости в одном заданном узле.
        """
        pass

    @abstractmethod
    def compute_aux_data_transient_field(self, time):
        """
        Вычисление вспомогательных данных зависимых от времени для моментов
        времени, в которое вычисляется все поле скорости.
        """
        pass

    def get_velocity_field(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Возвращает поле скорости на всей сетке."""
        return self._vel_field

    def get_pulsation_at_node(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Возвращает значение пульсаций в узле в заданные моменты времени."""
        return self._vel_puls

    @abstractmethod
    def compute_velocity_field(self, time):
        """Расчет поля скорости."""
        pass

    @abstractmethod
    def compute_pulsation_at_node(self):
        """Расчет пульсаций в заданном узле."""
        pass

    def get_puls_node(self):
        return self._i_puls, self._j_puls, self._k_puls

    def set_puls_node(self, i: int, j: int, k: int):
        self._i_puls = i
        self._j_puls = j
        self._k_puls = k

    @property
    def time_arr(self):
        return self._time_arr

    @time_arr.setter
    def time_arr(self, value: np.ndarray):
        self._ts = value[1] - value[0]
        self._num_ts = value.shape[0] - 1
        self._time_arr = value
