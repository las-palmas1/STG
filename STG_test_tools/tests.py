import unittest
from analysis_tools import Analyzer
from generators.storage import Lund, Smirnov, OriginalSEM, Davidson
from generators.abstract import Block, BCType
import numpy as np
import STG.common


class LundTest(unittest.TestCase):
    def setUp(self):
        n = 150
        # mesh = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), np.linspace(0, 1, n))
        mesh = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), [0, 1 / (n - 1)])
        self.block = Block(
            shape=(n, n, 2),
            mesh=(mesh[1], mesh[0], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        self.generator = Lund(
            block=self.block,
            u_av=(0., 0., 0.),
            re_uu=1.,
            re_vv=1.,
            re_ww=1.,
            re_uv=0.,
            re_uw=0.,
            re_vw=0.,
            time_arr=np.array([0])
        )
        self.analyzer = Analyzer(self.generator)

    def test_plot_2d_velocity_field(self):
        self.analyzer.plot_2d_velocity_field(vmin=-2.5, vmax=2.5, grid=False)

    def test_plot_velocity_history(self):
        self.analyzer.plot_velocity_history(0, 0, 0, 0.001, 1000)

    def test_plot_moments(self):
        self.analyzer.plot_moments(0, 0, 0, 0.001, 10000)

    def test_plot_divergence_field_2d(self):
        self.analyzer.plot_divergence_field_2d(grid=False)

    def test_plot_two_point_space_correlation(self):
        self.analyzer.plot_two_point_space_correlation(
            i0=0, j0=0, k0=0, ts=0.001, num_ts=1000, di=1, dj=1, dk=0, num=49
        )

    def test_plot_two_point_time_correlation(self):
        self.analyzer.plot_two_point_time_correlation(
            i=0, j=0, k=0, t0=0, t1=1.0, num_dt_av=500, num_dt=500
        )

    def test_plot_spectrum_2d(self):
        self.analyzer.plot_spectrum_2d(num_pnt=200)


class SmirnovTest(unittest.TestCase):
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
    def setUp(self):
        n = 55
        size = 20
        # mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), np.linspace(0, size, n))
        mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), [0, size / (n - 1)])
        self.block = Block(
            shape=(n, n, 2),
            mesh=(mesh[1], mesh[0], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        STG.common.init_rand()
        self.generator = Smirnov(
            block=self.block,
            u_av=(0., 0., 0.),
            ls_i=1.,
            ts_i=1,
            re_uu=1.,
            re_vv=1.,
            re_ww=1.,
            re_uv=0.,
            re_uw=0.,
            re_vw=0.,
            time_arr=np.array([0, 1]),
            mode_num=1000
        )
        self.analyzer = Analyzer(self.generator)

    def tearDown(self):
        self.analyzer.generator.free_data()

    def test_plot_2d_velocity_field(self):
        self.analyzer.plot_2d_velocity_field(vmin=-2.5, vmax=2.5, grid=False)

    def test_plot_velocity_history(self):
        self.analyzer.plot_velocity_history(0, 0, 0, 0.1, 1000)

    def test_plot_moments(self):
        self.analyzer.plot_moments(0, 0, 0, 0.01, 80000, ylim=(-0.5, 1.5))

    def test_plot_divergence_field_2d(self):
        self.analyzer.plot_divergence_field_2d(vmin=-1.5, vmax=1.5, grid=False)

    def test_plot_two_point_space_correlation(self):
        self.analyzer.plot_two_point_space_correlation(
            i0=0, j0=0, k0=0, ts=0.025, num_ts=4000, di=1, dj=1, dk=0, num=49
        )

    def test_plot_two_point_time_correlation(self):
        self.analyzer.plot_two_point_time_correlation(
            i=0, j=0, k=0, t0=0, t1=100, num_dt_av=4000, num_dt=200
        )

    def test_plot_spectrum_2d(self):
        self.analyzer.plot_spectrum_2d(num_pnt=200)

    def test_plot_spectrum_3d(self):
        self.analyzer.plot_spectrum_3d(num_pnt=200)


class DavidsonTest(unittest.TestCase):
    def tearDown(self):
        self.analyzer.generator.free_data()

    def setUp(self):
        n = 50
        size = 0.1 * n
        # mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), np.linspace(0, size, n))
        mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), [0, size / (n - 1)])
        self.block = Block(
            shape=(n, n, 2),
            mesh=(mesh[1], mesh[0], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        self.visc = 1.5e-5
        self.re_uu = 1.
        self.re_vv = 1.
        self.re_ww = 1.
        self.re_uv = 0
        self.re_uw = 0.
        self.re_vw = 0.
        self.ls_i = 2
        self.ts_i = 0.01
        self.dissip_rate = 0.09 ** 0.75 * (0.5 * (self.re_uu + self.re_vv + self.re_ww)) ** 1.5 / self.ls_i
        STG.common.init_rand()
        self.generator = Davidson(
            block=self.block,
            u_av=(0., 0., 0.),
            ls_i=self.ls_i,
            ts_i=self.ts_i,
            num_modes=1000,
            visc=self.visc,
            dissip_rate=self.dissip_rate,
            re_uu=self.re_uu,
            re_vv=self.re_vv,
            re_ww=self.re_ww,
            re_uv=self.re_uv,
            re_uw=self.re_uw,
            re_vw=self.re_vw,
            time_arr=np.array([0, 0.1]),
        )
        self.analyzer = Analyzer(self.generator)

    def test_plot_2d_velocity_field(self):
        self.analyzer.plot_2d_velocity_field(vmin=-2.5, vmax=2.5, grid=False)

    def test_plot_velocity_history(self):
        self.analyzer.plot_velocity_history(0, 0, 0, 0.001, 1000)

    def test_plot_moments(self):
        self.analyzer.plot_moments(0, 0, 0, 0.001, 5000, ylim=(-0.5, 1.5))

    def test_plot_divergence_field_2d(self):
        self.analyzer.plot_divergence_field_2d(vmin=-1.5, vmax=1.5, grid=False)

    def test_plot_two_point_space_correlation(self):
        self.analyzer.plot_two_point_space_correlation(
            i0=0, j0=0, k0=0, ts=0.002, num_ts=6000, di=1, dj=1, dk=0, num=49
        )

    def test_plot_two_point_time_correlation(self):
        self.analyzer.plot_two_point_time_correlation(
            i=0, j=0, k=0, t0=0, t1=1, num_dt_av=4000, num_dt=200
        )

    def test_plot_spectrum_2d(self):
        self.analyzer.plot_spectrum_2d(num_pnt=200)

    def test_plot_spectrum_3d(self):
        self.analyzer.plot_spectrum_3d(num_pnt=200)


class OriginalSEMTest(unittest.TestCase):
    def setUp(self):
        n = 70
        size = 6.28
        # mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), np.linspace(0, size, n))
        mesh = np.meshgrid(np.linspace(0, size, n), np.linspace(0, size, n), [0, size / (n - 1)])
        self.block = Block(
            shape=(n, n, 2),
            mesh=(mesh[1], mesh[0], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        self.generator = OriginalSEM(
            block=self.block,
            u_e=(10, 0, 0),
            re_uu=1.,
            re_vv=1.,
            re_ww=1.,
            re_uv=0.,
            re_uw=0.,
            re_vw=0.,
            time_arr=np.array([0]),
            sigma=0.5,
            eddies_num=1000
        )
        self.analyzer = Analyzer(self.generator)

    def test_plot_2d_velocity_field(self):
        self.analyzer.plot_2d_velocity_field(figsize=(7, 7), num_levels=20, vmin=-4, vmax=4, grid=False)

    def test_plot_velocity_history(self):
        self.analyzer.plot_velocity_history(10, 10, 0, 0.01, 1000)

    def test_plot_moments(self):
        self.analyzer.plot_moments(20, 20, 0, 0.005, 9000)

    def test_plot_divergence_field_2d(self):
        self.analyzer.plot_divergence_field_2d(vmin=-15, vmax=15, grid=False, num_levels=20)

    def test_plot_two_point_space_correlation(self):
        self.analyzer.plot_two_point_space_correlation(
            i0=0, j0=0, k0=0, ts=0.003, num_ts=6000, di=1, dj=1, dk=0, num=25
        )

    def test_plot_two_point_time_correlation(self):
        self.analyzer.plot_two_point_time_correlation(
            i=0, j=0, k=0, t0=0, t1=5, num_dt_av=1000, num_dt=1000)

    def test_plot_spectrum_2d(self):
        self.analyzer.plot_spectrum_2d(num_pnt=200)

    def test_plot_spectrum_3d(self):
        self.analyzer.plot_spectrum_3d(num_pnt=200)


