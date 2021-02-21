import pytest
from analysis_tools import Analyzer
from generators.storage import Lund, Smirnov, OriginalSEM, Davidson, Spectral
from generators.abstract import Block, BCType
import numpy as np
import STG.common
from STG.davidson import get_davidson_spectrum
import config
import config_plots
import os


def setup_module():
    if not os.path.exists(config_plots.plots_dir):
        os.makedirs(config_plots.plots_dir)


class Tester2D:
    def __init__(self, generator_type, **params):
        self.params = params
        mesh = STG.common.get_mesh(
            np.linspace(0, config.length, config.num_nodes),
            np.linspace(0, config.length, config.num_nodes),
            np.array([0, config.length / (config.num_nodes - 1)])
        )
        self.block = Block(
            shape=(config.num_nodes, config.num_nodes, 2),
            mesh=(mesh[0], mesh[1], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        if generator_type == Lund:
            self.method = 'Lund'
            self.generator = Lund(
                block=self.block,
                u_av=config.u_av,
                re_uu=config.re_uu,
                re_vv=config.re_vv,
                re_ww=config.re_ww,
                re_uv=config.re_uv,
                re_uw=config.re_uw,
                re_vw=config.re_vw,
                time_arr=np.array([0, 0.01 * config.length / max(config.u_av)])
            )
        elif generator_type == Smirnov:
            self.method = 'Smirnov'
            self.generator = Smirnov(
                block=self.block,
                u_av=config.u_av,
                re_uu=config.re_uu,
                re_vv=config.re_vv,
                re_ww=config.re_ww,
                re_uv=config.re_uv,
                re_uw=config.re_uw,
                re_vw=config.re_vw,
                time_arr=np.array([0, 0.01 * config.length / max(config.u_av)]),
                ls_i=config.ls_i,
                ts_i=config.ts_i,
                num_modes=params['num_modes']
            )
        elif generator_type == Davidson:
            self.method = 'Davidson'
            dissip_rate = 0.09 ** 0.75 * (0.5 * (config.re_uu + config.re_vv + config.re_ww)) ** 1.5 / config.ls_i
            visc = 1.5e-5
            self.generator = Davidson(
                block=self.block,
                u_av=config.u_av,
                re_uu=config.re_uu,
                re_vv=config.re_vv,
                re_ww=config.re_ww,
                re_uv=config.re_uv,
                re_uw=config.re_uw,
                re_vw=config.re_vw,
                time_arr=np.array([0, 0.01 * config.length / max(config.u_av)]),
                ls_i=config.ls_i,
                ts_i=config.ts_i,
                num_modes=params['num_modes'],
                dissip_rate=dissip_rate,
                visc=visc
            )
        elif generator_type == OriginalSEM:
            self.method = 'SEM'
            self.generator = OriginalSEM(
                block=self.block,
                u_e=config.u_av,
                re_uu=config.re_uu,
                re_vv=config.re_vv,
                re_ww=config.re_ww,
                re_uv=config.re_uv,
                re_uw=config.re_uw,
                re_vw=config.re_vw,
                time_arr=np.array([0, 0.01 * config.length / max(config.u_av)]),
                ls_ux=config.ls_ux, ls_uy=config.ls_uy, ls_uz=config.ls_uz,
                ls_vx=config.ls_vx, ls_vy=config.ls_vy, ls_vz=config.ls_vz,
                ls_wx=config.ls_wx, ls_wy=config.ls_wy, ls_wz=config.ls_wz,
                eddies_num=params['eddies_num']
            )
        elif generator_type == Spectral:
            dissip_rate = 0.09 ** 0.75 * (0.5 * (config.re_uu + config.re_vv + config.re_ww)) ** 1.5 / config.ls_i
            visc = 1.5e-5
            delta_min = STG.common.get_min_mesh_step(
                np.linspace(0, config.length, config.num_nodes),
                np.linspace(0, config.length, config.num_nodes),
                np.array([0, config.length / (config.num_nodes - 1)])
            )
            k_arr, energy = get_davidson_spectrum(
                delta_min, params['num_modes'], config.re_uu, config.re_vv, config.re_ww,
                config.ls_i, dissip_rate, visc
            )
            self.method = 'Spectral'
            self.generator = Spectral(
                block=self.block,
                u_av=config.u_av,
                re_uu=config.re_uu,
                re_vv=config.re_vv,
                re_ww=config.re_ww,
                ts_i=config.ts_i,
                num_modes=params['num_modes'],
                time_arr=np.array([0, 0.01 * config.length / max(config.u_av)]),
                k_arr=k_arr,
                energy=energy
            )
        else:
            self.method = ''
            self.generator = None
        STG.common.init_rand()
        self.analyzer = Analyzer(self.generator)

    @classmethod
    def dict_to_str(cls, d: dict):
        res = ''
        for key in d.keys():
            res += '-%s-%s' % (key, d[key])
        return res

    def plot_2d_velocity_field(self):
        self.analyzer.plot_2d_velocity_field(
            vmin=config_plots.vmin_vel_2d, vmax=config_plots.vmax_vel_2d,
            num_levels=config_plots.num_levels_vel_2d, grid=False, figsize=config_plots.figsize_contour,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.vel_2d_fname_base, method=self.method,
                                                params=self.dict_to_str(self.params))
            ),
            show=False, title_fsize=config_plots.title_fsize,
            axes=config_plots.axes_contour
        )

    def plot_velocity_history(self):
        num_tlvl = int(config.t_end_vel / config.t_step + 1)
        self.analyzer.plot_velocity_history(
            0, 0, 0, ts=config.t_step, num_tlvl=num_tlvl, ylim=config_plots.ylim_vel,
            figsize=config_plots.figsize_plot,
            label_fsize=config_plots.label_fsize, title_fsize=config_plots.title_fsize,
            ticks_fsize=config_plots.ticks_fsize, axes=config_plots.axes_plot,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.vel_hist_fname_base, method=self.method,
                                                params=self.dict_to_str(self.params))
            ),
            show=False
        )

    def plot_moments(self):
        num_tlvl = int(config.t_av / config.t_step + 1)
        self.analyzer.plot_moments(
            0, 0, 0, ts=config.t_step, num_tlvl=num_tlvl, figsize=config_plots.figsize_plot,
            ylim=config_plots.ylim_sec_mom, label_fsize=config_plots.label_fsize,
            ticks_fsize=config_plots.ticks_fsize, legend_fsize=config_plots.legend_fsize,
            axes=config_plots.axes_plot,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.moments_fname_base, method=self.method,
                                                params=self.dict_to_str(self.params))
            ),
            show=False
        )

    def plot_two_point_space_correlation(self, dir: str = 'x'):
        num_tlvl = int(config.t_av / config.t_step + 1)
        if dir == 'x':
            di = 1
            dj = 0
        elif dir == 'y':
            di = 0
            dj = 1
        else:
            di = 1
            dj = 1
        self.analyzer.plot_two_point_space_correlation(
            i0=0, j0=0, k0=0, ts=config.t_step, num_tlvl=num_tlvl, di=di, dj=dj, dk=0, num=config.num_nodes-1,
            figsize=config_plots.figsize_plot,
            ylim=config_plots.ylim_space_cor, label_fsize=config_plots.label_fsize,
            ticks_fsize=config_plots.ticks_fsize, legend_fsize=config_plots.legend_fsize,
            axes=config_plots.axes_plot,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.space_cor_fname_base, method=self.method,
                                                params='%s-dir-%s' % (self.dict_to_str(self.params), dir))
            ),
            show=False
        )

    def plot_two_point_time_correlation(self):
        t1 = 0.5 * config.t_av
        num_dt_av = int(t1 / config.t_step + 1)
        self.analyzer.plot_two_point_time_correlation(
            i=0, j=0, k=0, t0=0, t1=t1, num_dt_av=num_dt_av, num_dt=config.num_dt_acor,
            figsize=config_plots.figsize_plot,
            ylim=config_plots.ylim_auto_cor, label_fsize=config_plots.label_fsize,
            ticks_fsize=config_plots.ticks_fsize, legend_fsize=config_plots.legend_fsize,
            axes=config_plots.axes_plot,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.autocor_fname_base, method=self.method,
                                                params=self.dict_to_str(self.params))
            ),
            show=False
        )

    def plot_spectrum_2d(self):
        self.analyzer.plot_spectrum_2d(
            num_pnt=200, figsize=config_plots.figsize_spectrum, ylim=config_plots.ylim_spectrum,
            xlim=config_plots.xlim_spectrum, label_fsize=config_plots.label_fsize,
            ticks_fsize=config_plots.ticks_fsize, legend_fsize=config_plots.legend_fsize,
            axes=config_plots.axes_spectrum,
            fname=os.path.join(
                config_plots.plots_dir,
                config_plots.fname_templ.format(base=config_plots.spec_2d_fname_base, method=self.method,
                                                params=self.dict_to_str(self.params))
            ),
            show=False
        )


class TestLund2D:
    def setup(self):
        self.tester = Tester2D(Lund)

    def test_plot_2d_velocity_field(self):
        self.tester.plot_2d_velocity_field()

    def test_plot_velocity_history(self):
        self.tester.plot_velocity_history()

    def test_plot_moments(self):
        self.tester.plot_moments()

    @pytest.mark.parametrize(argnames='dir', argvalues=['x', 'y'])
    def test_plot_two_point_space_correlation(self, dir):
        self.tester.plot_two_point_space_correlation(dir)

    def test_plot_two_point_time_correlation(self):
        self.tester.plot_two_point_space_correlation()

    def test_plot_spectrum_2d(self):
        self.tester.plot_spectrum_2d()


class TestSmirnov2D:
    """
    При единичных временном и линейном масштабе метод Смирнова тождественнен методу Крайхмана. Линейный масштаб
    ни на что не влияет, зачем он нужен - неясно. Двумерный и трехмерный спектры сильно зависит от шага сетки,
    размера области, а также временного масштаба. Полного совпадения с заданным спектром  у двумерного спектра не
    выходит. При некоторых соотношении значений шага, размера области и временного масштаба удается добиться того,
    что к заданному спектру близка правая часть получаемого спектра в плоть до экстремума.
    Однако для трехмерного спектра добиться полного совпадения с заданным возможно.
    Дивергенция уменьшается с уменьшением шага. В пределе вроде бы получается нулевой.
    С графиками истории скорости в точки и графиками вторых моментов вроде бы все нормально. Вышли такими же,
    как и в статье, в которой изложен данный метод.
    """

    def teardown(self):
        self.tester.analyzer.generator.free_data()

    def setup(self):
        self.tester = Tester2D(Smirnov, num_modes=500)

    def test_plot_2d_velocity_field(self):
        self.tester.plot_2d_velocity_field()

    def test_plot_velocity_history(self):
        self.tester.plot_velocity_history()

    def test_plot_moments(self):
        self.tester.plot_moments()

    @pytest.mark.parametrize(argnames='dir', argvalues=['x', 'y'])
    def test_plot_two_point_space_correlation(self, dir):
        self.tester.plot_two_point_space_correlation(dir)

    def test_plot_two_point_time_correlation(self):
        self.tester.plot_two_point_time_correlation()

    def test_plot_spectrum_2d(self):
        self.tester.plot_spectrum_2d()


class TestDavidson2D:
    def teardown(self):
        self.tester.analyzer.generator.free_data()

    def setup(self):
        self.tester = Tester2D(Davidson, num_modes=200)

    def test_plot_2d_velocity_field(self):
        self.tester.plot_2d_velocity_field()

    def test_plot_velocity_history(self):
        self.tester.plot_velocity_history()

    def test_plot_moments(self):
        self.tester.plot_moments()

    @pytest.mark.parametrize(argnames='dir', argvalues=['x', 'y'])
    def test_plot_two_point_space_correlation(self, dir):
        self.tester.plot_two_point_space_correlation(dir)

    def test_plot_two_point_time_correlation(self):
        self.tester.plot_two_point_time_correlation()

    def test_plot_spectrum_2d(self):
        self.tester.plot_spectrum_2d()


class TestOriginalSEM2D:
    def teardown(self):
        self.tester.analyzer.generator.free_data()

    def setup(self):
        self.tester = Tester2D(OriginalSEM, eddies_num=300)

    def test_plot_2d_velocity_field(self):
        self.tester.plot_2d_velocity_field()

    def test_plot_velocity_history(self):
        self.tester.plot_velocity_history()

    def test_plot_moments(self):
        self.tester.plot_moments()

    @pytest.mark.parametrize(argnames='dir', argvalues=['x', 'y'])
    def test_plot_two_point_space_correlation(self, dir):
        self.tester.plot_two_point_space_correlation(dir)

    def test_plot_two_point_time_correlation(self):
        self.tester.plot_two_point_time_correlation()

    def test_plot_spectrum_2d(self):
        self.tester.plot_spectrum_2d()


class TestSpectral2D:
    def teardown(self):
        self.tester.analyzer.generator.free_data()

    def setup(self):
        self.tester = Tester2D(Spectral, num_modes=500)

    def test_plot_2d_velocity_field(self):
        self.tester.plot_2d_velocity_field()

    def test_plot_velocity_history(self):
        self.tester.plot_velocity_history()

    def test_plot_moments(self):
        self.tester.plot_moments()

    @pytest.mark.parametrize(argnames='dir', argvalues=['x', 'y'])
    def test_plot_two_point_space_correlation(self, dir):
        self.tester.plot_two_point_space_correlation(dir)

    def test_plot_two_point_time_correlation(self):
        self.tester.plot_two_point_time_correlation()

    def test_plot_spectrum_2d(self):
        self.tester.plot_spectrum_2d()
