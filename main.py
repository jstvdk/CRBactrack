import JF2012
from scipy.integrate import odeint
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm


def calc_new_coord(E, Z, l, b):
    q = -4.80320425e-10  # esu / antiproton charge
    Bx = 4.6e-6  # G
    lm = 1.  # kpc
    beta = 1.
    E = E  # eV
    Z = Z
    l = l * u.degree.to(u.rad)
    b = b * u.degree.to(u.rad)
    C1 = Z * q * (lm * 3.08567758e21) * Bx / (E * 1.60217733e-12)
    lc = 0.1
    Cb = (-Z * q / (E * 1.602e-12)) ** 2 * 2 / 9 * (lc / 3.24e-22)
    var0 = [beta * np.cos(l) * np.cos(b), beta * np.sin(l) * np.cos(b), beta * np.sin(b), -8.5, 0., 0., 0.]
    t = np.linspace(0, 200, 2000)

    def equation(var, t):
        bx = var[0]
        by = var[1]
        bz = var[2]
        x = var[3]
        y = var[4]
        z = var[5]
        dbxdt = C1 * (by * JF2012.Ballz(x, y, z) - bz * JF2012.Bally(x, y, z)) / Bx
        dbydt = -C1 * (bx * JF2012.Ballz(x, y, z) - bz * JF2012.Ballx(x, y, z)) / Bx
        dbzdt = C1 * (bx * JF2012.Bally(x, y, z) - by * JF2012.Ballx(x, y, z)) / Bx
        dxdt = bx
        dydt = by
        dzdt = bz
        dtheta2dt = Cb / 3.24e-22 * JF2012.BtotR(x, y, z) ** 2
        return [dbxdt, dbydt, dbzdt, dxdt, dydt, dzdt, dtheta2dt]

    result = odeint(equation, var0, t)

    beta_x_t = result[:, 0]
    beta_y_t = result[:, 1]
    beta_z_t = result[:, 2]
    x_t = result[:, 3]
    y_t = result[:, 4]
    z_t = result[:, 5]
    theta2_t = result[:, 6]


    def select_gal_prop(res):
        beta_x_t_gal = []
        beta_y_t_gal = []
        beta_z_t_gal = []
        x_t_gal = []
        y_t_gal = []
        z_t_gal = []
        theta2_t_gal = []
        for i in range(len(t)):
            if np.sqrt(res[:, 3][i] ** 2 + res[:, 4][i] ** 2 + res[:, 5][i] ** 2) < 20.:
                beta_x_t_gal.append(res[:, 0][i])
                beta_y_t_gal.append(res[:, 1][i])
                beta_z_t_gal.append(res[:, 2][i])
                x_t_gal.append(res[:, 3][i])
                y_t_gal.append(res[:, 4][i])
                z_t_gal.append(res[:, 5][i])
                theta2_t_gal.append(res[:, 6][i])
        coord = SkyCoord(x=beta_x_t_gal[-1], y=beta_y_t_gal[-1], z=beta_z_t_gal[-1],
                     unit='kpc', representation_type='cartesian')
        theta = np.sqrt(theta2_t_gal[-1]) * 360 / (2 * np.pi)
        return [coord, theta]
    new_lon = select_gal_prop(result)[0].spherical.lon.degree
    new_lat = select_gal_prop(result)[0].spherical.lat.degree
    theta = select_gal_prop(result)[1]
    return [E, new_lon, new_lat, beta_x_t, beta_y_t, beta_z_t, x_t, y_t, z_t, theta]


def plotSPtot():
    x = np.linspace(-20, 20, 800)
    y = np.linspace(-20, 20, 800)
    Bplot = np.array([[(JF2012.BSptot(i, j, 0.1) + JF2012.Bring(i, j, 0.1)) * 1e6 for i in x] for j in y])
    plt.contourf(x, y, Bplot, levels=500, cmap=cm.bwr)
    plt.colorbar()
    plt.clim(-4, 4)
    plt.show()
    return


def plot_trajectory(E, Z, l, b):
    fig, ax = plt.subplots()
    img = plt.imread('spiral_field.png')
    ax.set_ylim(-20, 20)
    ax.set_xlim(-20, 20)
    ax.imshow(img, extent=[-20, 20, -20, 20])
    ax.plot(calc_new_coord(E, Z, l, b)[6], calc_new_coord(E, Z, l, b)[7], 'g-')
    plt.grid()
    plt.show()


from Events import events_1020_no_calib, events_1020_ta1052_pao0948, events_1020_ta095_pao087


def write_new_events(events, z, filename):
    file = open('z{0}_{1}.txt'.format(z, filename), 'w+')
    for i in range(len(events)):
        file.write(str(events[i][0]) + ' ' +
                   str(calc_new_coord(events[i][0] * 1e18, z, events[i][1], events[i][2])[1]) + ' ' +
                   str(calc_new_coord(events[i][0] * 1e18, z, events[i][1], events[i][2])[2]) + ' ' +
                   str(calc_new_coord(events[i][0] * 1e18, z, events[i][1], events[i][2])[9]) + '\n')
    file.close()
    return 0.


write_new_events(events_1020_no_calib, z=6, filename='events_1020_no_calib_theta')
write_new_events(events_1020_ta1052_pao0948, z=6, filename='events_1020_ta1052_pao0948_theta')
write_new_events(events_1020_ta095_pao087, z=6, filename='events_1020_ta095_pao087_theta')

#plot_trajectory(1e20, 26, 35, -4)
# plotSPtot()