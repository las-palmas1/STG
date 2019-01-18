import unittest
from STG import common
from STG import smirnov
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
        mesh = np.meshgrid(np.linspace(1, 5, 5), np.linspace(6, 10, 5), np.linspace(11, 15, 5))
        self.re_uu = 1.
        self.re_vv = 2.
        self.re_ww = 3.
        self.init_data = common.get_init_data(
            mesh=(mesh[1], mesh[0], mesh[2]),
            re_uu=self.re_uu, re_vv=self.re_vv, re_ww=self.re_ww,
            re_uv=0., re_uw=0., re_vw=0., l_t=1., tau_t=1.
        )
        self.num_modes = 100
        self.time_arr = np.linspace(1, 5, 5)

    def test_compute_smirnov_data(self):
        data = smirnov.SmirnovData()
        smirnov.compute_Smirnov_data(self.init_data, self.num_modes, self.time_arr, data)
        self.assertEqual(self.re_uu, data.c1[1])
        self.assertEqual(self.re_vv, data.c2[2])
        self.assertEqual(self.re_ww, data.c3[0])
        self.assertNotEqual(data.p1[0], data.p1[2])
        self.assertNotEqual(data.p1[0], data.p2[0])
        self.assertNotEqual(data.p1[0], data.p3[0])
        self.assertNotEqual(data.p1[0], data.q1[0])
        smirnov.free_Smirnov_data(data)

    def test_compute_Smirnov_field_node(self):
        data = smirnov.SmirnovData()
        out_data = smirnov.OutDataNode()
        smirnov.compute_Smirnov_data(self.init_data, self.num_modes, self.time_arr, data)
        smirnov.compute_Smirnov_field_node(self.init_data, data, out_data, 0, 0, 0)
        u, v, w = smirnov.extract_pulsations_node(out_data)
        self.assertNotEqual(u[0], u[1])
        self.assertNotEqual(u[0], v[1])
        self.assertNotEqual(w[0], v[1])
        self.assertNotEqual(u[0], 0)
        self.assertIsNotNone(u[0])
        smirnov.free_Smirnov_data(data)
        smirnov.free_out_data_node(out_data)
