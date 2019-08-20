import unittest
from STG import common
from STG import smirnov
from STG import davidson
import matplotlib.pyplot as plt
import numpy as np


class TestProdDistSamples(unittest.TestCase):
    def setUp(self):
        common.init_rand()

    def test_uniform(self):
        x1 = -5
        x2 = 5
        num = 2000
        samp = common.get_uniform(x1, x2, num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(x1, x2, 20))
        plt.grid()
        plt.show()

    def test_normal(self):
        mu = 1
        sigma = 1
        num = 2000
        samp = common.get_normal(mu, sigma, num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(-2, 4, 20))
        plt.grid()
        plt.show()

    def test_trigon(self):
        num = 5000
        samp = common.get_trigon(num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(0, np.pi, 20))
        plt.grid()
        plt.show()


class TestSmirnov(unittest.TestCase):
    def setUp(self):
        mesh = np.meshgrid(np.linspace(1, 5, 2), np.linspace(6, 10, 7), np.linspace(11, 15, 1))
        self.re_uu = 1.
        self.re_vv = 2.
        self.re_ww = 3.
        self.init_data = common.get_init_data(
            mesh=(mesh[1], mesh[0], mesh[2]),
            re_uu=self.re_uu, re_vv=self.re_vv, re_ww=self.re_ww,
            re_uv=0., re_uw=0., re_vw=0.,
            ls_i=1., ls_ux=0., ls_uy=0., ls_uz=0.,
            ls_vx=0., ls_vy=0., ls_vz=0.,
            ls_wx=0., ls_wy=0., ls_wz=0.,
            ts_i=1., ts_u=0., ts_v=0., ts_w=0.
        )
        self.num_modes = 100
        self.time_arr = np.linspace(1, 5, 5)

    def test_compute_smirnov_data(self):
        data = smirnov.STG_SmirnovData()
        smirnov.compute_smirnov_data(self.init_data, self.num_modes, data)
        self.assertAlmostEqual(self.re_uu**0.5, data.c1[1], places=5)
        self.assertAlmostEqual(self.re_vv**0.5, data.c2[2], places=5)
        self.assertAlmostEqual(self.re_ww**0.5, data.c3[0], places=5)
        self.assertNotEqual(data.p1[0], data.p1[2])
        self.assertNotEqual(data.p1[0], data.p2[0])
        self.assertNotEqual(data.p1[0], data.p3[0])
        self.assertNotEqual(data.p1[0], data.q1[0])
        smirnov.free_smirnov_data(data)

    def test_compute_Smirnov_field_node(self):
        data = smirnov.STG_SmirnovData()
        node_hist = smirnov.STG_VelNodeHist()
        smirnov.compute_smirnov_data(self.init_data, self.num_modes, data)
        num_ts = len(self.time_arr - 1)
        ts = self.time_arr[1] - self.time_arr[0]
        smirnov.compute_smirnov_node_hist(self.init_data, data, ts, num_ts, node_hist, 0, 0, 0)
        u, v, w = common.extract_pulsations_from_node_hist(node_hist)
        self.assertNotEqual(u[0], u[1])
        self.assertNotEqual(u[0], v[1])
        self.assertNotEqual(w[0], v[1])
        self.assertNotEqual(u[0], 0)
        self.assertIsNotNone(u[0])
        smirnov.free_smirnov_data(data)
        smirnov.free_node_hist(node_hist)


class TestDavidson(unittest.TestCase):
    def setUp(self):
        mesh = np.meshgrid(np.linspace(0.1, 0.5, 2), np.linspace(0.6, 0.10, 7), np.linspace(0.11, 0.15, 2))
        self.re_uu = 1.
        self.re_vv = 2.
        self.re_ww = 3.
        self.ls_i = 0.02
        self.ts_i = 1.
        self.init_data = common.get_init_data(
            mesh=(mesh[1], mesh[0], mesh[2]),
            re_uu=self.re_uu, re_vv=self.re_vv, re_ww=self.re_ww,
            re_uv=0., re_uw=0., re_vw=0.,
            ls_i=self.ls_i, ls_ux=0., ls_uy=0., ls_uz=0.,
            ls_vx=0., ls_vy=0., ls_vz=0.,
            ls_wx=0., ls_wy=0., ls_wz=0.,
            ts_i=self.ts_i, ts_u=0., ts_v=0., ts_w=0.
        )
        self.num_modes = 100
        self.time_arr = np.linspace(1, 5, 5)
        self.visc = 1.8e-5
        self.dissip_rate = 0.09**0.75 * (0.5 * (self.re_uu + self.re_vv + self.re_ww))**1.5 / self.ls_i

    def test_compute_davidson_stat_data(self):
        ts = self.time_arr[1] - self.time_arr[0]
        stat_data = davidson.STG_DavidsonData_Stationary()
        davidson.compute_davidson_stat_data(self.init_data, self.num_modes, self.dissip_rate,
                                            self.visc, ts, stat_data)
        trace = self.re_uu + self.re_vv + self.re_ww
        self.assertAlmostEqual((self.re_uu*3 / trace)**0.5, stat_data.c1[1], places=5)
        self.assertAlmostEqual((self.re_vv*3 / trace)**0.5, stat_data.c2[2], places=5)
        self.assertAlmostEqual((self.re_ww*3 / trace)**0.5, stat_data.c3[0], places=5)
        davidson.free_davidson_stat_data(stat_data)

    def test_compute_davidson_node_hist(self):
        ts = self.time_arr[1] - self.time_arr[0]
        stat_data = davidson.STG_DavidsonData_Stationary()
        davidson.compute_davidson_stat_data(self.init_data, self.num_modes, self.dissip_rate,
                                            self.visc, ts, stat_data)
        trans_data = davidson.STG_DavidsonData_Transient()
        davidson.alloc_davidson_trans_data(self.init_data, self.num_modes, trans_data)
        node_hist = common.STG_VelNodeHist()
        davidson.compute_davidson_node_hist(
            self.init_data, stat_data, ts, len(self.time_arr) - 1, trans_data, node_hist, 0, 0, 0)
        u, v, w = common.extract_pulsations_from_node_hist(node_hist)
        self.assertNotEqual(u[0], u[1])
        self.assertNotEqual(u[0], v[1])
        self.assertNotEqual(w[0], v[1])
        self.assertNotEqual(u[0], 0)
        self.assertIsNotNone(u[0])
        common.free_node_hist(node_hist)
        davidson.free_davidson_stat_data(stat_data)
        davidson.free_davidson_trans_data(trans_data)



