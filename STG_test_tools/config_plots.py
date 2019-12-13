
# путь к папке с картинками
plots_dir = 'plots'

vel_2d_fname_base = 'vel_2d'
vel_hist_fname_base = 'vel_hist'
moments_fname_base = 'moments'
space_cor_fname_base = 'space_cor'
autocor_fname_base = 'autocor'
spec_2d_fname_base = 'spec_2d'
spec_3d_fname_base = 'spec_3d'

fname_templ = '{base:s}-{method:s}{params:s}'

vmin_vel_2d = -25
vmax_vel_2d = 25
num_levels_vel_2d = 30

# настроки осей для всех графиков, кроме профиля сокрости и спектра
axes_plot = (0.11, 0.14, 0.86, 0.84)
# настройки осей для профиля сокрости
axes_contour = (0.05, 0.05, 0.9, 0.9)
# настройки осей для спектра
axes_spectrum = (0.15, 0.11, 0.82, 0.86)

# размеры шрифтов
legend_fsize = 14
ticks_fsize = 14
label_fsize = 16
title_fsize = 18

# пределы по вертикальной оси для различных графиков
ylim_vel = (-25, 25)
ylim_sec_mom = (-10, 140)
ylim_space_cor = (-0.4, 1.1)
ylim_auto_cor = (-0.4, 1.1)
ylim_spectrum = (1e-4, 1e1)

# пределы по горизонтальной оси
xlim_spectrum = (1e-1, 1e3)

# размеры картинок
figsize_contour = (6, 6)
figsize_plot = (6, 4)
figsize_spectrum = (6, 6)
