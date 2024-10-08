import JF2012
from scipy.integrate import odeint
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm

# Const

q = -4.80320425e-10  # esu / antiproton charge
Bx = 4.6e-6  # G
lm = 1.  # kpc
E = 1e20  # eV
Z = 26.
l = 35 * u.degree.to(u.rad)
b = 10 * u.degree.to(u.rad)
C1 = Z * q * (lm * 3.08567758e21) * Bx / (E * 1.60217733e-12)
beta = 1.
lc = 0.1
Cb = (-Z * q / (E * 1.602e-12)) ** 2 * 2 / 9 * (lc / 3.24e-22)


def calc_new_coord(data):
    return 0.

def plotSPtot():
    x = np.linspace(-20, 20, 800)
    y = np.linspace(-20, 20, 800)
    Bplot = np.array([[(JF2012.BSptot(i, j, 0.1) + JF2012.Bring(i, j, 0.1)) * 1e6 for i in x] for j in y])
    plt.contourf(x, y, Bplot, levels=500, cmap=cm.bwr)
    plt.colorbar()
    plt.clim(-4, 4)
    plt.show()
    return


def equation(u, t):
    bx = u[0]
    by = u[1]
    bz = u[2]
    x = u[3]
    y = u[4]
    z = u[5]
    #theta2 = u[6]
    dbxdt = C1 * (by * JF2012.Ballz(x, y, z) - bz * JF2012.Bally(x, y, z)) / Bx
    dbydt = -C1 * (bx * JF2012.Ballz(x, y, z) - bz * JF2012.Ballx(x, y, z)) / Bx
    dbzdt = C1 * (bx * JF2012.Bally(x, y, z) - by * JF2012.Ballx(x, y, z)) / Bx
    dxdt = bx
    dydt = by
    dzdt = bz
    dtheta2dt = Cb / 3.24e-22 * JF2012.BtotR(x, y, z) ** 2
    return [dbxdt, dbydt, dbzdt, dxdt, dydt, dzdt, dtheta2dt]


u0 = [beta * np.cos(l) * np.cos(b), beta * np.sin(l) * np.cos(b), beta * np.sin(b), -8.5, 0., 0., 0.]

t = np.linspace(0, 20, 2000)
print(t)

u = odeint(equation, u0, t)
beta_x_t = u[:, 0]
beta_y_t = u[:, 1]
beta_z_t = u[:, 2]
x_t = u[:, 3]
y_t = u[:, 4]
z_t = u[:, 5]
theta2_t = u[:, 6]


beta_x_t_gal = []
beta_y_t_gal = []
beta_z_t_gal = []
x_t_gal = []
y_t_gal = []
z_t_gal = []
theta2_t_gal = []

print(np.sqrt(theta2_t) * 360 / (2 * np.pi))

def select_gal_prop(u):
    for i in range(len(t)):
        if np.sqrt(u[:, 3][i] ** 2 + u[:, 4][i] ** 2 + u[:, 5][i] ** 2) < 20.:
            beta_x_t_gal.append(u[:, 0][i])
            beta_y_t_gal.append(u[:, 1][i])
            beta_z_t_gal.append(u[:, 2][i])
            x_t_gal.append(u[:, 3][i])
            y_t_gal.append(u[:, 4][i])
            z_t_gal.append(u[:, 5][i])
            theta2_t_gal.append(u[:, 6][i])
            # c_lon.append(SkyCoord(x=beta_x_t_gal[i], y=beta_y_t_gal[i], z=beta_z_t_gal[i],
            #                       unit='kpc', representation_type='cartesian').spherical.lon.degree)
            # c_lat.append(SkyCoord(x=beta_x_t_gal[i], y=beta_y_t_gal[i], z=beta_z_t_gal[i],
            #                       unit='kpc', representation_type='cartesian').spherical.lat.degree)
        #     print(np.sqrt(u[:, 3][i] ** 2 + u[:, 4][i] ** 2 + u[:, 5][i] ** 2))
        #     print('bingo')
        # else:
        #     print(np.sqrt(u[:, 3][i] ** 2 + u[:, 4][i] ** 2 + u[:, 5][i] ** 2))
        #     print('false')

new_lon = []
new_lat = []


select_gal_prop(u)
print(np.sqrt(theta2_t_gal[-1]) * 360 / (2 * np.pi))

c = SkyCoord(x=beta_x_t_gal[-1], y=beta_y_t_gal[-1], z=beta_z_t_gal[-1],
             unit='kpc', representation_type='cartesian')
c2 = SkyCoord(x=np.cos(l) * np.cos(b), y=np.sin(l) * np.cos(b), z=np.sin(b), unit='kpc', representation_type='cartesian')
new_lon.append(c.spherical.lon.degree)
new_lat.append(c.spherical.lat.degree)
print(c.spherical)
print(c2.spherical)
print(t[-1], beta_x_t_gal[-1], beta_y_t_gal[-1], beta_z_t_gal[-1])
print(np.sqrt(x_t_gal[-1] ** 2 + y_t_gal[-1] ** 2 + z_t_gal[-1] ** 2))
print(new_lon)
print(new_lat)

fig, ax = plt.subplots()
img = plt.imread('spiral_field.png')
ax.set_ylim(-20, 20)
ax.set_xlim(-20, 20)
ax.imshow(img, extent=[-20, 20, -20, 20])
ax.plot(x_t, y_t, 'g-')
plt.grid()

plt.show()


# plotSPtot()