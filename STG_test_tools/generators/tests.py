from generators.storage import OriginalSEM, Block, BCType
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as ax3
import matplotlib.animation as anim
import unittest
from STG.sem import STG_SEMData_Stationary, STG_SEMData_Transient
from STG.common import get_mesh, init_rand


def get_limits(stat_data: STG_SEMData_Stationary):
    x_min = stat_data.vol_lims.x_min
    x_max = stat_data.vol_lims.x_max
    y_min = stat_data.vol_lims.y_min
    y_max = stat_data.vol_lims.y_max
    z_min = stat_data.vol_lims.z_min
    z_max = stat_data.vol_lims.z_max
    return x_min, x_max, y_min, y_max, z_min, z_max


def get_eddies_pos(trans_data: STG_SEMData_Transient):
    num_ts = trans_data.num_ts
    num_eddies = trans_data.num_eddies
    x_e = np.zeros([num_ts + 1, num_eddies])
    y_e = np.zeros([num_ts + 1, num_eddies])
    z_e = np.zeros([num_ts + 1, num_eddies])

    for i_ts in range(num_ts + 1):
        for i_e in range(num_eddies):
            x_e[i_ts, i_e] = trans_data.eddies_pos[i_ts * num_eddies + i_e].x
            y_e[i_ts, i_e] = trans_data.eddies_pos[i_ts * num_eddies + i_e].y
            z_e[i_ts, i_e] = trans_data.eddies_pos[i_ts * num_eddies + i_e].z
    return x_e, y_e, z_e


class TestSEM(unittest.TestCase):
    def setUp(self):
        n = 150
        size = 1
        mesh = get_mesh(np.linspace(0, size, n), np.linspace(0, size, n), np.linspace(0, size, n))
        self.block = Block(
            shape=(n, n, n),
            mesh=(mesh[0], mesh[1], mesh[2]),
            bc=[(BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall), (BCType.NotWall, BCType.NotWall)]
        )
        l_i = 0.05

        init_rand()
        self.sem = OriginalSEM(
            block=self.block,
            u_e=(1, 0, 0),
            re_uu=1.,
            re_vv=1.,
            re_ww=1.,
            re_uv=0.,
            re_uw=0.,
            re_vw=0.,
            ls_ux=l_i, ls_uy=l_i, ls_uz=l_i,
            ls_vx=l_i, ls_vy=l_i, ls_vz=l_i,
            ls_wx=l_i, ls_wy=l_i, ls_wz=l_i,
            eddies_num=100,
            time_arr=np.linspace(0, 1, 4)
        )
        limits = get_limits(self.sem.c_stat_data)
        self.x_min = limits[0]
        self.x_max = limits[1]
        self.y_min = limits[2]
        self.y_max = limits[3]
        self.z_min = limits[4]
        self.z_max = limits[5]
        self.x_e, self.y_e, self.z_e = get_eddies_pos(self.sem.c_trans_data)
        self.num_time_levels = self.sem.c_trans_data.num_ts + 1

    def tearDown(self):
        self.sem.free_data()

    def test_plot_eddies_centers_movement(self):
        fig = plt.figure(figsize=(9, 7))
        ax = ax3.Axes3D(fig)

        ax.plot(xs=[self.x_min, self.x_min], ys=[self.y_min, self.y_min],
                zs=[self.z_min, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_min], ys=[self.y_min, self.y_max],
                zs=[self.z_min, self.z_min], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_max], ys=[self.y_min, self.y_min],
                zs=[self.z_min, self.z_min], c='red', lw=2)
        ax.plot(xs=[self.x_max, self.x_max], ys=[self.y_min, self.y_min],
                zs=[self.z_min, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_max, self.x_max], ys=[self.y_min, self.y_max],
                zs=[self.z_min, self.z_min], c='red', lw=2)
        ax.plot(xs=[self.x_max, self.x_max], ys=[self.y_max, self.y_max],
                zs=[self.z_min, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_max], ys=[self.y_max, self.y_max],
                zs=[self.z_min, self.z_min], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_min], ys=[self.y_max, self.y_max],
                zs=[self.z_min, self.z_max], c='red', lw=2)

        ax.plot(xs=[self.x_min, self.x_max], ys=[self.y_min, self.y_min],
                zs=[self.z_max, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_min], ys=[self.y_min, self.y_max],
                zs=[self.z_max, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_max, self.x_max], ys=[self.y_min, self.y_max],
                zs=[self.z_max, self.z_max], c='red', lw=2)
        ax.plot(xs=[self.x_min, self.x_max], ys=[self.y_max, self.y_max],
                zs=[self.z_max, self.z_max], c='red', lw=2)
        line = ax.plot(self.x_e[0, :], self.y_e[0, :], self.z_e[0, :], ls='', marker='o')[0]

        def update(frame):
            line.set_data([self.x_e[frame, :], self.y_e[frame, :]])
            line.set_3d_properties(self.z_e[frame, :])

        ani = anim.FuncAnimation(fig, func=update, frames=self.num_time_levels,
                                 interval=400)
        plt.show()

