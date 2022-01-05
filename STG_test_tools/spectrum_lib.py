import numpy as np


def _get_energy(m_grid: np.ndarray, energy_arr: np.ndarray, m_mag):
    """Вычисляет величину энергии в шаровом слое единичной толщины в простанстве m"""
    energy_arr_filt = energy_arr[(m_grid > m_mag - 0.5) * (m_grid < m_mag + 0.5)]
    energy = energy_arr_filt.sum()
    return energy


def get_spectrum_1d(x: np.ndarray, y: np.ndarray, band_width: int = 5, band_growth='const', **opts):
    """
    Для определения спектральной плотности интвервал волновых чисел (k_min, k_max), соответствующий дискр-му
    преобразованию Фурье F{y} разделяется на полосы. На каждой такой полосе находится средняя спектральная
    плотность путем суммирования энергии мод и отнесении этой величины к размеру полосы.

    Для определения величин полос есть два способа. При первом способе (band_growth='const') все полосы одинаковы
    по ширине и включают себя band_width элементов . При втором способе (band_growth='exp') ширина полос
    увеличивается с ростом частоты в степенном масштабе. Такой рост обеспечивается следующим образом. Вводится
    новая величина k_l(k), равная k_l(k) = k**exp, где 0 < exp < 1. Интервал (k_l(k_min), k_l(k_max)) разбивается
    на равные полосы таким образом, что первая полоса в пространстве k включает в себя band_width элементов.
    """
    if len(x) != len(y):
        msg = "Can't compute spectrum: x and y have different sizes: '{}' and '{}'.".format(len(x), len(y))
        raise ValueError(msg)
    if band_width <= 0:
        msg = "band_width must be positive"
        raise ValueError(msg)

    band_growth_values = ('exp', 'const')
    if band_growth not in band_growth_values:
        msg = "band_growth has value '{}', but it was expected to be one " \
              "of '{}'".format(band_growth, band_growth_values)
        raise ValueError(msg)

    if band_growth == 'exp':
        exp = opts.setdefault('exp', 0.5)
        if not 0 < exp < 1:
            msg = "exp is equal to '{}', but it was expected to be greater than 0 and less than 1.".format(exp)
            raise ValueError(msg)
    else:
        exp = None

    num = len(x)
    dx = x[1] - x[0]
    length = dx * num
    # шаг волнового числа (циклической частоты в случае, если x - время)
    dk = 2 * np.pi / length
    # шаг частоты в Гц (в случае, если x - время)
    df = 1 / length
    if num % 2 == 0:
        k_indexes_p = np.linspace(0, num / 2 - 1, num // 2, dtype=np.int)
        k_indexes_n = np.linspace(-num / 2, -1, num // 2, dtype=np.int)
    else:
        k_indexes_p = np.linspace(0, (num - 1) / 2, (num - 1) // 2 + 1, dtype=np.int)
        k_indexes_n = np.linspace(-(num - 1) / 2, -1, (num - 1) // 2, dtype=np.int)

    # Расположение индексов в порядке, используемом в numpy и scipy (сначала неотрицательные,
    # а затем отр-ые частоты в порядке роста частот)
    k_indexes = np.append(k_indexes_p, k_indexes_n)

    k = dk * k_indexes
    y_image = np.fft.fft(y) / num

    # оставляем только положительные частоты
    k_filter = np.where(k > 0)[0]
    num_cut = len(k_filter)
    k_cut = k[k_filter]
    y_image_cut = y_image[k_filter]
    # энергия мод
    e_modes_cut = np.real(y_image_cut * np.conjugate(y_image_cut))

    if band_growth == 'const':
        band_num = num_cut // band_width
        # границы полос
        k_bands = np.linspace(k_cut.min(), band_width * dk * band_num, band_num + 1)
    elif band_growth == 'exp':
        k_log = k_cut**exp
        dk_band_log = k_cut[band_width]**exp - k_cut[0]**exp
        band_num = int((np.max(k_log) - np.min(k_log)) // dk_band_log)
        k_log_bands = np.linspace(np.min(k_log), dk_band_log * band_num, band_num + 1)
        k_bands = k_log_bands**(1/exp)
    else:
        band_num = num_cut // band_width
        # границы полос
        k_bands = np.linspace(k_cut.min(), band_width * dk * band_num, band_num + 1)

    # определение спектральной плотности энергии
    k_spec = k_bands[0: band_num]
    f_spec = k_spec / dk * df
    e_spec = np.zeros(band_num)
    for i in range(band_num):
        e_modes_band = np.sum(e_modes_cut[np.where((k_cut >= k_bands[i]) * (k_cut < k_bands[i + 1]))[0]])
        e_spec[i] = 2 * e_modes_band / (k_bands[i + 1] - k_bands[i])
    return k_spec, f_spec, e_spec


def get_spectrum_2d(x: np.ndarray, vel: np.ndarray, num_pnt=100):
    """Вычисление двухмерного энергетического спектра для компоненты скорости на квадратной сетке."""
    vel_f = np.fft.fftn(vel)
    vel_f = np.fft.fftshift(vel_f)
    length = x.max() - x.min()
    m1, m2 = np.meshgrid(
        np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
        np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
        indexing='ij'
    )
    m_grid = np.sqrt(m1**2 + m2**2)
    energy = 0.5 * (np.abs(vel_f) / (x.shape[0]**2)) ** 2

    m_mag = np.linspace(0, m_grid.max(), num_pnt)
    e_k_mag = np.zeros(num_pnt)
    for i in range(num_pnt):
        e_k_mag[i] = _get_energy(m_grid, energy, m_mag[i]) * length / (2 * np.pi)
    k_mag = m_mag * 2 * np.pi / length
    return k_mag, e_k_mag


def get_spectrum_3d(x: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray, num_pnt=100):
    """Вычисление двухмерного энергетического спектра для компоненты скорости на квадратной сетке."""
    u_f = np.fft.fftn(u)
    v_f = np.fft.fftn(v)
    w_f = np.fft.fftn(w)
    u_f = np.fft.fftshift(u_f)
    v_f = np.fft.fftshift(v_f)
    w_f = np.fft.fftshift(w_f)
    length = x.max() - x.min()
    m1, m2, m3 = np.meshgrid(
        np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
        np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
        np.linspace(-x.shape[0] / 2, x.shape[0] / 2 - 1, x.shape[0]),
        indexing='ij'
    )
    m_grid = np.sqrt(m1 ** 2 + m2 ** 2 + m3 ** 2)
    energy = 0.5 * (np.abs(u_f)**2 + np.abs(v_f)**2 + np.abs(w_f)**2) * (1 / x.shape[0])**6
    m_mag = np.linspace(0, m_grid.max(), num_pnt)
    e_k_mag = np.zeros(num_pnt)
    for i in range(num_pnt):
        e_k_mag[i] = _get_energy(m_grid, energy, m_mag[i]) * length / (2 * np.pi)
    k_mag = m_mag * 2 * np.pi / length
    return k_mag, e_k_mag