
STG_lib_name = 'libstg'

# TODO: сделать возможность генерировать неоднородную сетку, задавая распределения вдоль каждой из осей
#  добавить возможность задавать гетерогенность параметров вдоль заданной оси (например, ввести в качестве
#  параметра ось (x или y) и задавать в нормализованном виде распределения или считывать их из файлов)

# Число узлов
num_nodes = 50
# Размер области
length = 0.25
u_av = (480, 0, 0)

# Шаг по времени
t_step = 1e-5
# Интервал для проб скорости
t_end_vel = 1 * 0.1
# Интервал осреднения для расчета статистики
t_av = 1 * 0.1
# Число шагов между концами интервала, на котором расчитывается автокорреляции
num_dt_acor = 50

# Рейнольдсовы напряжения
re_uu = 24**2
re_vv = 24**2
re_ww = 24**2
re_uv = 0
re_uw = 0
re_vw = 0

# Масштабы
k_ls = 0.07 / length
ls_ux = k_ls * length
ls_uy = k_ls * length
ls_uz = k_ls * length
ls_vx = k_ls * length
ls_vy = k_ls * length
ls_vz = k_ls * length
ls_wx = k_ls * length
ls_wy = k_ls * length
ls_wz = k_ls * length
ls_i = (ls_ux + ls_uy + ls_vx + ls_vy + ls_wx + ls_wy) / 9

k_ts = 1 / t_step * k_ls * length / u_av[0]
ts_u = k_ts * t_step
ts_v = k_ts * t_step
ts_w = k_ts * t_step
ts_i = (ts_u + ts_v + ts_w) / 3
