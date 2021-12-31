from generators.abstract import Generator
import matplotlib.pyplot as plt
from typing import Tuple
import numpy as np
import enum
from scipy.optimize import fsolve


class Spec1DMode(enum.Enum):
    F_MODE = 0
    K_MODE = 1


class Analyzer:
    def __init__(self, generator: Generator):
        self.generator = generator
        self.generator.compute_velocity_field(num_ts=0)

    def plot_2d_velocity_field(
            self, figsize=(7, 7), num_levels=20, vmin=-3.5, vmax=3.5, grid=True, title_fsize=18, title=True,
            axes=(0.05, 0.05, 0.9, 0.9), fname=None, show=False
    ):
        x = self.generator.block.mesh[0][:, :, 0]
        y = self.generator.block.mesh[1][:, :, 0]
        u3d, v3d, w3d = self.generator.get_velocity_field()
        u = u3d[:, :, 0]
        v = v3d[:, :, 0]
        w = w3d[:, :, 0]

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.contourf(x, y, u, num_levels, cmap='rainbow', vmin=vmin, vmax=vmax)
        if title:
            plt.title('U', fontsize=title_fsize, fontweight='bold')
        plt.colorbar()
        plt.xticks(x[:, 0], [])
        plt.yticks(y[0, :], [])
        plt.grid(grid)
        if fname:
            plt.savefig(fname + '-U')
        if show:
            plt.show()

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.contourf(x, y, v, num_levels, cmap='rainbow', vmin=vmin, vmax=vmax)
        if title:
            plt.title('V', fontsize=title_fsize, fontweight='bold')
        plt.colorbar()
        plt.xticks(x[:, 0], [])
        plt.yticks(y[0, :], [])
        plt.grid(grid)
        if fname:
            plt.savefig(fname + '-V')
        if show:
            plt.show()

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.contourf(x, y, w, num_levels, cmap='rainbow', vmin=vmin, vmax=vmax)
        if title:
            plt.title('W', fontsize=title_fsize, fontweight='bold')
        plt.colorbar()
        plt.xticks(x[:, 0], [])
        plt.yticks(y[0, :], [])
        plt.grid(grid)
        if fname:
            plt.savefig(fname + '-W')
        if show:
            plt.show()

    def save_velocity_field_tec(self, fname):
        pass

    def plot_velocity_history(
            self, i: int, j: int, k: int, ts, num_tlvl: int, figsize=(7, 7), ylim=(-3, 3),
            label_fsize=16, ticks_fsize=12, title_fsize=18, title=True, axes=(0.07, 0.07, 0.9, 0.87), fname=None,
            show=False
    ):
        t_arr = np.arange(0, ts * num_tlvl, ts)
        self.generator.set_puls_node(i, j, k)
        self.generator.time_arr = t_arr
        self.generator.compute_pulsation_at_node()
        u, v, w = self.generator.get_pulsation_at_node()

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(t_arr, u, color='red', lw=1)
        plt.grid()

        if ylim:
            plt.ylim(*ylim)
        plt.xlim(xmin=0, xmax=ts * num_tlvl)
        if title:
            plt.title('U', fontsize=title_fsize, fontweight='bold')
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('t, с', fontsize=label_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        if fname:
            plt.savefig(fname + '-U')
        if show:
            plt.show()

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(t_arr, v, color='red', lw=1)
        plt.grid()

        if ylim:
            plt.ylim(*ylim)
        plt.xlim(xmin=0, xmax=ts * num_tlvl)
        if title:
            plt.title('V', fontsize=title_fsize, fontweight='bold')
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('t, с', fontsize=label_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        if fname:
            plt.savefig(fname + '-V')
        if show:
            plt.show()

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(t_arr, w, color='red', lw=1)
        plt.grid()

        if ylim:
            plt.ylim(*ylim)
        plt.xlim(xmin=0, xmax=ts * num_tlvl)
        if title:
            plt.title('W', fontsize=title_fsize, fontweight='bold')
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('t, с', fontsize=label_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        if fname:
            plt.savefig(fname + '-W')
        if show:
            plt.show()

    @classmethod
    def _get_average_arr(cls, value: np.ndarray):
        res = []
        sum = 0
        for n, i in enumerate(value):
            sum += i
            res.append(sum / (n + 1))
        return np.array(res)

    @classmethod
    def _get_average(cls, value: np.ndarray):
        return value.sum() / len(value)

    def plot_moments(
            self, i: int, j: int, k: int, ts: float, num_tlvl: int, figsize=(7, 7), ylim=(-0.5, 1),
            legend_fsize=14, ticks_fsize=12, label_fsize=14, axes=(0.07, 0.07, 0.9, 0.87), fname=None, show=False
    ):
        t_arr = np.arange(0, ts * num_tlvl, ts)
        self.generator.set_puls_node(i, j, k)
        self.generator.time_arr = t_arr
        self.generator.compute_pulsation_at_node()
        u, v, w = self.generator.get_pulsation_at_node()
        uu_av = self._get_average_arr(u * u)
        vv_av = self._get_average_arr(v * v)
        ww_av = self._get_average_arr(w * w)
        uv_av = self._get_average_arr(u * v)
        uw_av = self._get_average_arr(u * w)
        vw_av = self._get_average_arr(v * w)

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(t_arr, uu_av, lw=1.5, color='red', label=r'$<v_x^2>$')
        plt.plot(t_arr, vv_av, lw=1.5, color='blue', label=r'$<v_y^2>$')
        plt.plot(t_arr, ww_av, lw=1.5, color='green', label=r'$<v_z^2>$')
        plt.plot(t_arr, uv_av, lw=1.5, color='red', ls=':', label=r'$<v_x v_y>$')
        plt.plot(t_arr, uw_av, lw=1.5, color='blue', ls=':', label=r'$<v_x v_z>$')
        plt.plot(t_arr, vw_av, lw=1.5, color='green', ls=':', label=r'$<v_y v_z>$')

        plt.legend(fontsize=legend_fsize)
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('t, с', fontsize=label_fsize, fontweight='bold')
        plt.xlim(0, num_tlvl * ts)
        if ylim:
            plt.ylim(*ylim)
        plt.grid()
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()

    def plot_divergence_field_2d(
            self, figzie=(7, 7), num_levels=20, vmin=-300, vmax=300, grid=True, title_fsize=18,
            fname=None, show=False
    ):
        x = self.generator.block.mesh[0][:, :, 0]
        y = self.generator.block.mesh[1][:, :, 0]
        vel = self.generator.get_velocity_field()
        div = self.generator.get_divergence(vel, self.generator.block.mesh, self.generator.block.shape)
        div_2d = div[:, :, 0]

        plt.figure(figsize=figzie)
        plt.axes([0.05, 0.05, 0.9, 0.9])
        plt.contourf(x, y, div_2d, num_levels, cmap='rainbow', vmin=vmin, vmax=vmax)
        plt.title(r'$div(U)$', fontsize=title_fsize, fontweight='bold')
        plt.colorbar()
        plt.xticks(x[:, 0], [])
        plt.yticks(y[0, :], [])
        plt.grid(grid)
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()

    def save_divergence_field_tec(self, fname):
        pass

    def plot_two_point_space_correlation(
            self, i0: int, j0: int, k0: int, ts: float, num_tlvl: int,
            di: int = 1, dj: int = 1, dk: int = 1, num: int = 20,
            figsize=(7, 7), label_fsize=14, ticks_fsize=12, legend_fsize=14, ylim=(-1.1, 1.1),
            axes=(0.07, 0.07, 0.9, 0.87), fname=None, show=False
    ):
        t_arr = np.arange(0, ts * num_tlvl, ts)
        r = np.zeros(num)
        cor_uu = np.zeros(num)
        cor_vv = np.zeros(num)
        cor_ww = np.zeros(num)
        self.generator.set_puls_node(i0, j0, k0)
        self.generator.time_arr = t_arr
        self.generator.compute_pulsation_at_node()
        u0, v0, w0 = self.generator.get_pulsation_at_node()
        u0u0_av = self._get_average(u0 * u0)
        v0v0_av = self._get_average(v0 * v0)
        w0w0_av = self._get_average(w0 * w0)

        for n in range(0, num):
            i = i0 + di * n
            j = j0 + dj * n
            k = k0 + dk * n
            self.generator.set_puls_node(i, j, k)
            self.generator.compute_pulsation_at_node()
            u, v, w = self.generator.get_pulsation_at_node()
            r[n] = np.sqrt((self.generator.block.mesh[0][i, j, k] - self.generator.block.mesh[0][i0, j0, k0])**2 +
                           (self.generator.block.mesh[1][i, j, k] - self.generator.block.mesh[1][i0, j0, k0])**2 +
                           (self.generator.block.mesh[2][i, j, k] - self.generator.block.mesh[2][i0, j0, k0])**2)
            uu_av = self._get_average(u * u)
            vv_av = self._get_average(v * v)
            ww_av = self._get_average(w * w)
            u0u_av = self._get_average(u0 * u)
            v0v_av = self._get_average(v0 * v)
            w0w_av = self._get_average(w0 * w)
            cor_uu[n] = u0u_av / (np.sqrt(u0u0_av) * np.sqrt(uu_av))
            cor_vv[n] = v0v_av / (np.sqrt(v0v0_av) * np.sqrt(vv_av))
            cor_ww[n] = w0w_av / (np.sqrt(w0w0_av) * np.sqrt(ww_av))

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(r, cor_uu, color='red', lw=1.5, label=r'$R_{xx}^r$')
        plt.plot(r, cor_vv, color='blue', lw=1.5, label=r'$R_{yy}^r$')
        plt.plot(r, cor_ww, color='green', lw=1.5, label=r'$R_{zz}^r$')

        plt.grid()
        plt.xlim(xmin=0, xmax=r.max())
        if ylim:
            plt.ylim(*ylim)
        plt.xlabel('r, м', fontsize=label_fsize, fontweight='bold')
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        plt.legend(fontsize=legend_fsize)
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()

    def plot_two_point_time_correlation(
            self, i: int, j: int, k: int, t1: float, t0: float=0.,
            num_dt_av: int=200, num_dt: int=100, figsize=(6.5, 4.5),
            label_fsize=14, ticks_fsize=12, legend_fsize=14, ylim=(-1.1, 1.1),
            axes=(0.07, 0.07, 0.9, 0.87), fname=None, show=False
    ):
        """
        :param num_dt - число отрезков между моментами t1 и t2.
        :param num_dt_av - число отрезков между моментами t0 и t1 и оно же - половина числа
            отрезков ни интервале осреднения.

        [t1, t2] - интервал, на котором рассчитываются автокорреляции. Он разбит на num_dt отрезков.

        1. Момент t2 определяется как t1 + num_dt * (t1 - t0) / num_dt_av.
        2. Момент t3 определяется как t2 + (t1 - t0).
        3. Интервал [t2, t3] при этом разбивается на num_dt_av отрезков.
        4. Пульсации считаются во всех точках между моментами t0 и t3.
        5. Интервал осреднения T = t1 - t0, шаг dt = (t1 - t0) / num_dt_av.
        6. Осреднение проводится от момента t1 - T + i * num_dt_av до момента
            t1 + T + i * num_dt_av, где i = 0, 1, ... num_dt.
        """
        t2 = t1 + num_dt * (t1 - t0) / num_dt_av
        t3 = t2 + t1 - t0
        t_arr = np.linspace(t0, t3, num_dt_av * 2 + num_dt + 1)
        dt_arr = np.linspace(t1, t2, num_dt + 1) - t1
        cor_uu = np.zeros(num_dt + 1)
        cor_vv = np.zeros(num_dt + 1)
        cor_ww = np.zeros(num_dt + 1)
        self.generator.set_puls_node(i, j, k)
        self.generator.time_arr = t_arr
        self.generator.compute_pulsation_at_node()
        u, v, w = self.generator.get_pulsation_at_node()
        u0u0_av = self._get_average(u[0: 2 * num_dt_av + 1] * u[0: 2 * num_dt_av + 1])
        v0v0_av = self._get_average(v[0: 2 * num_dt_av + 1] * v[0: 2 * num_dt_av + 1])
        w0w0_av = self._get_average(w[0: 2 * num_dt_av + 1] * w[0: 2 * num_dt_av + 1])
        for n in range(num_dt + 1):
            uu_av = self._get_average(u[n: 2 * num_dt_av + n + 1] * u[n: 2 * num_dt_av + n + 1])
            vv_av = self._get_average(v[n: 2 * num_dt_av + n + 1] * v[n: 2 * num_dt_av + n + 1])
            ww_av = self._get_average(w[n: 2 * num_dt_av + n + 1] * w[n: 2 * num_dt_av + n + 1])
            u0u_av = self._get_average(u[0: 2 * num_dt_av + 1] * u[n: 2 * num_dt_av + n + 1])
            v0v_av = self._get_average(v[0: 2 * num_dt_av + 1] * v[n: 2 * num_dt_av + n + 1])
            w0w_av = self._get_average(w[0: 2 * num_dt_av + 1] * w[n: 2 * num_dt_av + n + 1])
            cor_uu[n] = u0u_av / (np.sqrt(u0u0_av) * np.sqrt(uu_av))
            cor_vv[n] = v0v_av / (np.sqrt(v0v0_av) * np.sqrt(vv_av))
            cor_ww[n] = w0w_av / (np.sqrt(w0w0_av) * np.sqrt(ww_av))

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(dt_arr, cor_uu, color='red', lw=1.5, label=r'$R_{xx}^t$')
        plt.plot(dt_arr, cor_vv, color='blue', lw=1.5, label=r'$R_{yy}^t$')
        plt.plot(dt_arr, cor_ww, color='green', lw=1.5, label=r'$R_{zz}^t$')

        plt.grid()
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlim(xmin=0, xmax=dt_arr.max())
        if ylim:
            plt.ylim(*ylim)
        plt.xlabel(r't, с', fontsize=label_fsize, fontweight='bold')
        plt.legend(fontsize=legend_fsize)
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()

    @classmethod
    def _get_energy(cls, m_grid: np.ndarray, energy_arr: np.ndarray, m_mag):
        """Вычисляет величину энергии в шаровом слое единичной толщины в простанстве m"""
        energy_arr_filt = energy_arr[(m_grid > m_mag - 0.5) * (m_grid < m_mag + 0.5)]
        energy = energy_arr_filt.sum()
        return energy

    @classmethod
    def _get_spectrum_1d(cls, x: np.ndarray, y: np.ndarray, band_width: int = 5, band_growth='const', **opts):
        """
        Для определения спектральной плотности интвервал волновых чисел (k_min, k_max), соответствующий дискр-му
        преобразованию Фурье F{y} разделяется на полосы. На каждой такой полосе находится средняя спектральная
        плотность путем суммирования энергии мод и отнесении этой величины к размеру полосы.

        Для определения величин полос есть два способа. При первом способе (band_growth='const') все полосы одинаковы
        по ширине и включают себя band_width элементов . При втором способе (band_growth='exp') ширина полос
        увеличивается с ростом частоты в степенном масштабе. Такой рост обеспечивается следующим образом. Вводится
        новая величина k_l(k), равная k_l(k) = k**exp, где 0 < exp < 1. Интервал (k_l(k_min), k_l(k_max)) разбивается
        на равные полосы таким образом, что первая полоса в пространстве k включает в себя band_width элементов.
        """
        if len(x) != len(y):
            msg = "Can't compute spectrum: x and y have different sizes: '{}' and '{}'.".format(len(x), len(y))
            raise ValueError(msg)
        if band_width <= 0:
            msg = "band_width must be positive"
            raise ValueError(msg)

        band_growth_values = ('exp', 'const')
        if band_growth not in band_growth_values:
            msg = "band_growth has value '{}', but it was expected to be one " \
                  "of '{}'".format(band_growth, band_growth_values)
            raise ValueError(msg)

        if band_growth == 'exp':
            exp = opts.setdefault('exp', 0.5)
            if not 0 < exp < 1:
                msg = "exp is equal to '{}', but it was expected to be greater than 0 and less than 1.".format(exp)
                raise ValueError(msg)
        else:
            exp = None

        num = len(x)
        dx = x[1] - x[0]
        length = dx * num
        # шаг волнового числа (циклической частоты в случае, если x - время)
        dk = 2 * np.pi / length
        # шаг частоты в Гц (в случае, если x - время)
        df = 1 / length
        if num % 2 == 0:
            k_indexes_p = np.linspace(0, num / 2 - 1, num // 2, dtype=np.int)
            k_indexes_n = np.linspace(-num / 2, -1, num // 2, dtype=np.int)
        else:
            k_indexes_p = np.linspace(0, (num - 1) / 2, (num - 1) // 2 + 1, dtype=np.int)
            k_indexes_n = np.linspace(-(num - 1) / 2, -1, (num - 1) // 2, dtype=np.int)

        # Расположение индексов в порядке, используемом в numpy и scipy (сначала неотрицательные,
        # а затем отр-ые частоты в порядке роста частот)
        k_indexes = np.append(k_indexes_p, k_indexes_n)

        k = dk * k_indexes
        y_image = np.fft.fft(y) / num

        # оставляем только положительные частоты
        k_filter = np.where(k > 0)[0]
        num_cut = len(k_filter)
        k_cut = k[k_filter]
        y_image_cut = y_image[k_filter]
        # энергия мод
        e_modes_cut = np.real(y_image_cut * np.conjugate(y_image_cut))

        if band_growth == 'const':
            band_num = num_cut // band_width
            # границы полос
            k_bands = np.linspace(k_cut.min(), band_width * dk * band_num, band_num + 1)
        elif band_growth == 'exp':
            k_log = k_cut**exp
            dk_band_log = k_cut[band_width]**exp - k_cut[0]**exp
            band_num = int((np.max(k_log) - np.min(k_log)) // dk_band_log)
            k_log_bands = np.linspace(np.min(k_log), dk_band_log * band_num, band_num + 1)
            k_bands = k_log_bands**(1/exp)
        else:
            band_num = num_cut // band_width
            # границы полос
            k_bands = np.linspace(k_cut.min(), band_width * dk * band_num, band_num + 1)

        # определение спектральной плотности энергии
        k_spec = k_bands[0: band_num]
        f_spec = k_spec / dk * df
        e_spec = np.zeros(band_num)
        for i in range(band_num):
            e_modes_band = np.sum(e_modes_cut[np.where((k_cut >= k_bands[i]) * (k_cut < k_bands[i + 1]))[0]])
            e_spec[i] = 2 * e_modes_band / (k_bands[i + 1] - k_bands[i])
        return k_spec, f_spec, e_spec

    @classmethod
    def _get_spectrum_2d(cls, x: np.ndarray, vel: np.ndarray, num_pnt=100):
        """Вычисление двухмерного энергетического спектра для компоненты скорости на квадратной сетке."""
        vel_f = np.fft.fftn(vel)
        vel_f = np.fft.fftshift(vel_f)
        length = x.max() - x.min()
        m1, m2 = np.meshgrid(
            np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
            np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
            indexing='ij'
        )
        m_grid = np.sqrt(m1**2 + m2**2)
        energy = 0.5 * (np.abs(vel_f) / (x.shape[0]**2)) ** 2

        m_mag = np.linspace(0, m_grid.max(), num_pnt)
        e_k_mag = np.zeros(num_pnt)
        for i in range(num_pnt):
            e_k_mag[i] = cls._get_energy(m_grid, energy, m_mag[i]) * length / (2 * np.pi)
        k_mag = m_mag * 2 * np.pi / length
        return k_mag, e_k_mag

    @classmethod
    def _get_spectrum_3d(cls, x: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray, num_pnt=100):
        """Вычисление двухмерного энергетического спектра для компоненты скорости на квадратной сетке."""
        u_f = np.fft.fftn(u)
        v_f = np.fft.fftn(v)
        w_f = np.fft.fftn(w)
        u_f = np.fft.fftshift(u_f)
        v_f = np.fft.fftshift(v_f)
        w_f = np.fft.fftshift(w_f)
        length = x.max() - x.min()
        m1, m2, m3 = np.meshgrid(
            np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
            np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
            np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
            indexing='ij'
        )
        m_grid = np.sqrt(m1 ** 2 + m2 ** 2 + m3 ** 2)
        energy = 0.5 * (np.abs(u_f)**2 + np.abs(v_f)**2 + np.abs(w_f)**2) * (1 / x.shape[0])**6
        m_mag = np.linspace(0, m_grid.max(), num_pnt)
        e_k_mag = np.zeros(num_pnt)
        for i in range(num_pnt):
            e_k_mag[i] = cls._get_energy(m_grid, energy, m_mag[i]) * length / (2 * np.pi)
        k_mag = m_mag * 2 * np.pi / length
        return k_mag, e_k_mag

    def compute_spectrum_1d(
            self, i: int, j: int, k: int, ts: float, num_tlvl: int,  band_width: int = 5,
            band_growth='exp',  **opts
    ):
        t_arr = np.arange(0, ts * num_tlvl, ts)
        self.generator.set_puls_node(i, j, k)
        self.generator.time_arr = t_arr
        self.generator.compute_pulsation_at_node()
        u, v, w = self.generator.get_pulsation_at_node()

        k_spec_u, f_spec_u, e_spec_u = self._get_spectrum_1d(t_arr, u, band_width, band_growth, **opts)
        k_spec_v, f_spec_v, e_spec_v = self._get_spectrum_1d(t_arr, v, band_width, band_growth, **opts)
        k_spec_w, f_spec_w, e_spec_w = self._get_spectrum_1d(t_arr, w, band_width, band_growth, **opts)
        return (k_spec_u, f_spec_u, e_spec_u), (k_spec_v, f_spec_v, e_spec_v), (k_spec_w, f_spec_w, e_spec_w),

    @classmethod
    def _get_53_line(cls, k_min: float, k_max: float, k0: float, e0: float, num=100):
        """
        Рассчитывает параметры линии e = c0 * k^(-5/3), проходящей через заданную точку и
        находящейся в заданном интервале волновых чисел.
        """
        c0 = fsolve(lambda c: e0 - c * k0**(-5 / 3), x0=np.array([1]))
        k = np.linspace(k_min, k_max, num)
        f = k / (2 * np.pi)
        e = c0 * k**(-5 / 3)
        return k, f, e

    @classmethod
    def plot_spectrum_1d(
            cls, u_spec: Tuple[np.ndarray, np.ndarray, np.ndarray],
            v_spec: Tuple[np.ndarray, np.ndarray, np.ndarray],
            w_spec: Tuple[np.ndarray, np.ndarray, np.ndarray],
            spec_mode: Spec1DMode = Spec1DMode.F_MODE,
            plot_all: bool = False,
            axes_rect=(0.15, 0.13, 0.83, 0.85),
            tick_fontsize=14,
            label_fontsize=16,
            legend_fontsize=18,
            ax_extend=0.08,
            fname=''
    ):
        figsize = (8, 6)
        lw = 1.5
        legend_loc = 'upper right'

        fig: plt.Figure = plt.Figure(figsize=figsize)
        ax: plt.Axes = fig.add_axes(axes_rect)

        if spec_mode == Spec1DMode.F_MODE:
            x_u = u_spec[1]
            x_v = u_spec[1]
            x_w = u_spec[1]
            xlabel = r'$f,\ Гц$'
        elif spec_mode == Spec1DMode.K_MODE:
            x_u = u_spec[0]
            x_v = u_spec[0]
            x_w = u_spec[0]
            xlabel = r'$\omega,\ 1/с$'
        else:
            x_u = u_spec[1]
            x_v = u_spec[1]
            x_w = u_spec[1]
            xlabel = r'$f,\ Гц$'

        if plot_all:
            x_min = min(np.min(x_u), np.min(x_v), np.min(x_w))
            x_max = min(np.max(x_u), np.max(x_v), np.max(x_w))
            y_min = max(np.min(u_spec[2]), np.min(v_spec[2]), np.min(w_spec[2]))
            y_max = max(np.max(u_spec[2]), np.max(v_spec[2]), np.max(w_spec[2]))
        else:
            x_min = np.min(x_u)
            x_max = np.max(x_u)
            y_min = np.min(u_spec[2])
            y_max = np.max(u_spec[2])

        dx_log = np.log10(x_max / x_min)
        x_min_ax_log = np.log10(x_min) - dx_log * ax_extend
        x_min_ax = 10**x_min_ax_log
        x_max_ax_log = np.log10(x_max) + dx_log * ax_extend
        x_max_ax = 10**x_max_ax_log

        dy_log = np.log10(y_max / y_min)
        y_min_ax_log = np.log10(y_min) - dy_log * ax_extend
        y_min_ax = 10**y_min_ax_log
        y_max_ax_log = np.log10(y_max) + dy_log * ax_extend
        y_max_ax = 10**y_max_ax_log

        if spec_mode == Spec1DMode.K_MODE:
            k_53, f_53, e_53 = cls._get_53_line(
                x_min_ax, x_max_ax, k0=0.5 * (x_min_ax + x_max_ax),
                e0=0.6 * (y_min_ax + y_max_ax)
            )
            x_53 = k_53
        elif spec_mode == Spec1DMode.F_MODE:
            k_53, f_53, e_53 = cls._get_53_line(
                x_min_ax * 2 * np.pi, x_max_ax * 2 * np.pi, k0=0.5 * (x_min_ax + x_max_ax) * 2 * np.pi,
                e0=0.6 * (y_min_ax + y_max_ax)
            )
            x_53 = f_53
        else:
            k_53, f_53, e_53 = cls._get_53_line(
                x_min_ax, x_max_ax, k0=0.5 * (x_min_ax + x_max_ax),
                e0=0.6 * (y_min_ax + y_max_ax)
            )
            x_53 = k_53

        ax.plot(x_u, u_spec[2], lw=lw, color='tab:red', label=r'$E_u$')
        ax.plot(x_53, e_53, lw=lw, color='green', label=r'$e = C k^{-5/3}$')
        if plot_all:
            ax.plot(x_v, v_spec[2], lw=lw, color='tab:blue', label=r'$E_v$')
            ax.plot(x_w, w_spec[2], lw=lw, color='tab:green', label=r'$E_w$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid()
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(r'$E,\ м^2/с$', fontsize=label_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        ax.legend(fontsize=legend_fontsize, loc=legend_loc)
        ax.set_xlim(x_min_ax, x_max_ax)
        ax.set_ylim(y_min_ax, y_max_ax)

        fig.savefig(fname)

    def plot_spectrum_2d(
            self, figsize=(7, 7), num_pnt=100, ylim=(1e-3, 1e1), xlim=None,
            label_fsize=14, ticks_fsize=12, legend_fsize=14,
            axes=(0.07, 0.07, 0.9, 0.87), fname=None, show=False
    ):
        i_num = self.generator.block.shape[0]
        j_num = self.generator.block.shape[1]
        if i_num != j_num:
            msg = "Mesh block has unequal numbers of nodes of its edges: {}, {}. " \
                  "2-D spectrum can not be computed.".format(i_num, j_num)
            raise Exception(msg)

        x = self.generator.block.mesh[0][:, 0, 0]
        y = self.generator.block.mesh[1][0, :, 0]
        dx = np.max(x) - np.min(x)
        dy = np.max(y) - np.min(y)
        if dx != dy:
            msg = "Block domain has edges of unequal length: {:.2e}, {:.2e}. " \
                  "2-D spectrum can not be computed.".format(dx, dy)
            raise Exception(msg)

        u3d, v3d, w3d = self.generator.get_velocity_field()
        u = u3d[:, :, 0]
        v = v3d[:, :, 0]
        w = w3d[:, :, 0]
        k_u, e_u = self._get_spectrum_2d(x, u, num_pnt)
        k_v, e_v = self._get_spectrum_2d(x, v, num_pnt)
        k_w, e_w = self._get_spectrum_2d(x, w, num_pnt)

        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(k_u, e_u, color='red', label=r'$E_u$', lw=1.5)
        plt.plot(k_v, e_v, color='blue', label=r'$E_v$', lw=1.5)
        plt.plot(k_w, e_w, color='green', label=r'$E_w$', lw=1.5)
        plt.plot(k_w, e_u + e_v + e_w, color='black', label=r'$E_\Sigma$', lw=2.5)

        if xlim:
            k = np.logspace(int(np.log10(xlim[0])), int(np.log10(xlim[1])), 500)
        else:
            k = np.logspace(-2, 3, 500)
        if (self.generator.get_desired_spectrum(k) == 0).all():
            pass
        else:
            plt.plot(k, self.generator.get_desired_spectrum(k), color='black', ls='--', lw=1.5, label='Заданный')

        if ylim:
            plt.ylim(*ylim)
        if xlim:
            plt.xlim(*xlim)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which='both')
        plt.legend(fontsize=legend_fsize)
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('k, 1/м', fontsize=label_fsize, fontweight='bold')
        plt.ylabel('E, м^3/с^2', fontsize=label_fsize, fontweight='bold')
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()

    def plot_spectrum_3d(
            self, figsize=(7, 7), num_pnt=100, ylim=(1e-3, 1e1), xlim=None,
            label_fsize=14, ticks_fsize=12, legend_fsize=14,
            axes=(0.07, 0.07, 0.9, 0.87), fname=None, show=False
    ):
        i_num = self.generator.block.shape[0]
        j_num = self.generator.block.shape[1]
        k_num = self.generator.block.shape[2]
        if i_num != j_num or i_num != k_num:
            msg = "Mesh block has unequal numbers of nodes of its edges: {}, {}, {}. " \
                  "3-D spectrum can not be computed.".format(i_num, j_num, k_num)
            raise Exception(msg)

        x = self.generator.block.mesh[0][:, 0, 0]
        y = self.generator.block.mesh[1][0, :, 0]
        z = self.generator.block.mesh[2][0, 0, :]
        dx = np.max(x) - np.min(x)
        dy = np.max(y) - np.min(y)
        dz = np.max(z) - np.min(z)
        if dx != dy or dx != dz:
            msg = "Block domain has edges of unequal length: {:.2e}, {:.2e}, {:.2e}. " \
                  "3-D spectrum can not be computed.".format(dx, dy, dz)
            raise Exception(msg)

        u, v, w = self.generator.get_velocity_field()
        k, e = self._get_spectrum_3d(x, u, v, w, num_pnt)
        plt.figure(figsize=figsize)
        if axes:
            plt.axes(axes)
        plt.plot(k, e, color='red', lw=2.5, label=r'$E_\Sigma$', )

        if xlim:
            k = np.logspace(int(np.log10(xlim[0])), int(np.log10(xlim[1])), 500)
        else:
            k = np.logspace(-2, 3, 500)
        if (self.generator.get_desired_spectrum(k) == 0).all():
            pass
        else:
            plt.plot(k, self.generator.get_desired_spectrum(k), color='black', ls='--', lw=1.5, label='Заданный')

        k_53 = np.logspace(-2, 2, 100)
        plt.plot(0.5 * k_53**(-5 / 3), k_53, color='black', ls=':', lw=1.5, label=r'$\sim k^{-\frac{5}{3}}$')

        if ylim:
            plt.ylim(*ylim)
        if xlim:
            plt.xlim(*xlim)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which='both')
        plt.legend(fontsize=legend_fsize)
        plt.xticks(fontsize=ticks_fsize, fontweight='bold')
        plt.yticks(fontsize=ticks_fsize, fontweight='bold')
        plt.xlabel('k, 1/м', fontsize=label_fsize, fontweight='bold')
        plt.ylabel('E, м^3/с^2', fontsize=label_fsize, fontweight='bold')
        if fname:
            plt.savefig(fname)
        if show:
            plt.show()















