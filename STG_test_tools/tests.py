import unittest
import SGT
import matplotlib.pyplot as plt
import numpy as np


class TestProdDistSamples(unittest.TestCase):

    def test_uniform(self):
        x1 = -5
        x2 = 5
        num = 2000
        samp = SGT.get_uniform(x1, x2, num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(x1, x2, 20))
        plt.grid()
        plt.show()

    def test_normal(self):
        mu = 1
        sigma = 1
        num = 2000
        samp = SGT.get_normal(mu, sigma, num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(-2, 4, 20))
        plt.grid()
        plt.show()

    def test_trigon(self):
        num = 5000
        samp = SGT.get_trigon(num)
        plt.figure(figsize=(7, 5))
        plt.hist(samp, bins=np.linspace(0, np.pi, 20))
        plt.grid()
        plt.show()
