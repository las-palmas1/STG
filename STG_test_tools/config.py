
STG_lib_name = 'libstg'

# TODO: сделать возможность генерировать неоднородную сетку, задавая распределения вдоль каждой из осей
#  добавить возможность задавать гетерогенность параметров вдоль заданной оси (например, ввести в качестве
#  параметра ось (x или y) и задавать в нормализованном виде распределения или считывать их из файлов)

# Число узлов
num_nodes_x = 51
num_nodes_y = 51
# Размер области
length_x = 2.2 * 0.1
# length_y = length_x / (num_nodes_x - 1) * (num_nodes_y - 1)
length_y = 3 * 0.1
u_av = (481.2, 0, 0)

# Масштабы
k_ls = 0.05 / length_x
# ls_ux = k_ls * length_x
# ls_uy = k_ls * length_x
# ls_uz = k_ls * length_x
# ls_vx = k_ls * length_x
# ls_vy = k_ls * length_x
# ls_vz = k_ls * length_x
# ls_wx = k_ls * length_x
# ls_wy = k_ls * length_x
# ls_wz = k_ls * length_x
l_t = 0.01
ls_ux = l_t
ls_uy = l_t
ls_uz = l_t
ls_vx = l_t
ls_vy = l_t
ls_vz = l_t
ls_wx = l_t
ls_wy = l_t
ls_wz = l_t
ls_i = (ls_ux + ls_uy + ls_vx + ls_vy + ls_wx + ls_wy) / 6

u_abs = (u_av[0]**2 + u_av[1]**2 + u_av[2]**2)**0.5

ts_u = ls_i / u_abs
ts_v = ls_i / u_abs
ts_w = ls_i / u_abs
ts_i = (ts_u + ts_v + ts_w) / 3


# Шаг по времени
t_step = ts_i / 10
# t_step = 1e-6
# t_step = 0.1 / 1000
# Интервал для проб скорости
# t_end_vel = 1000 * t_step
t_end_vel = 0.1 / 20
# Интервал осреднения для расчета статистики
# t_av = 5000 * t_step
t_av = 0.1
# Число шагов между концами интервала, на котором расчитывается автокорреляции
num_dt_acor = 50

# Рейнольдсовы напряжения
re_uu = 579.039
re_vv = 579.039
re_ww = 579.039
re_uv = 0
re_uw = 0
re_vw = 0
