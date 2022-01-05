
# путь к папке с картинками
plots_dir = 'plots-ONERA-It-0.05-Lt-0.01'

vel_2d_fname_base = 'vel_2d'
vel_hist_fname_base = 'vel_hist'
moments_fname_base = 'moments'
space_cor_fname_base = 'space_cor'
autocor_fname_base = 'autocor'
spec_1d_fname_base = 'spec_1d-dt_rel-10'
spec_2d_fname_base = 'spec_2d'
spec_3d_fname_base = 'spec_3d'

fname_templ = '{base:s}-{method:s}{params:s}'

vmin_vel_2d = -40
vmax_vel_2d = 40
num_levels_vel_2d = 41

spec_1d_band_width = 4
spec_1d_band_growth = 'exp'
spec_1d_exp = 0.6

# настроки осей для всех графиков, кроме профиля сокрости и спектра
axes_plot = (0.12, 0.14, 0.85, 0.84)
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
ylim_vel = (-100, 100)
ylim_sec_mom = (-1, 2)
ylim_space_cor = (-0.4, 1.1)
ylim_auto_cor = (-0.4, 1.1)
ylim_spectrum_23d = (1e-2, 1e1)
ylim_spectrum_1d = (1e-9, 1e0)

# пределы по горизонтальной оси
xlim_spectrum_23d = (1e-1, 1e3)
xlim_spectrum_1d = (1e0, 5e5)

# размеры картинок
figsize_contour = (6, 6)
figsize_plot = (6, 4)
figsize_spectrum = (6, 5)
