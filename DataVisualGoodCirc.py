from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Events import events_1020_no_calib, events_1020_ta095_pao087, events_1020_ta1052_pao0948
from Sources import AGNs, SBGs
import os
#os.environ["PROJ_LIB"] = "D:\\Anaconda3\\envs\\plot\\Library\\share"  # fix for basemap
from mpl_toolkits.basemap import Basemap


def read_events(filename):
    # Функція для читання з файлу нових положень КП та розмазування у випадковому маг. полі
    file = open(filename, 'r')
    l = []
    b = []
    theta = []
    for line in file.readlines():
        l.append(float(line.split()[1]))
        b.append(float(line.split()[2]))
        theta.append(float(line.split()[3]))
    return np.array([l, b, theta])


init_evs_l = np.array([events_1020_ta095_pao087[i][1] for i in range(len(events_1020_ta095_pao087))])  # Початкові полження КП: l
init_evs_b = np.array([events_1020_ta095_pao087[i][2] for i in range(len(events_1020_ta095_pao087))])  # Початкові полження КП: b

evs1 = read_events('z1_events_1020_ta095_pao087_theta.txt')  # Події і кути розмазування для z=1
evs6 = read_events('z6_events_1020_ta095_pao087_theta.txt')  # Події і кути розмазування для z=6
evs26 = read_events('z26_events_1020_ta095_pao087_theta.txt')  # Події і кути розмазування для z=26


# SOURCES start_
starburst = SBGs
AGNs = AGNs


def circle(alpha, l0, b0):
    alpha = np.radians(alpha)
    l0 = np.radians(l0)
    b0 = np.radians(-(90 - b0))
    t = np.linspace(0, 2 * np.pi, 400)
    x = np.sin(alpha) * np.cos(b0) * np.cos(l0) * np.cos(t) + np.sin(alpha) * np.sin(l0) * np.sin(t) - np.cos(alpha) * np.sin(b0) * np.cos(l0)
    y = -np.sin(alpha) * np.cos(b0) * np.sin(l0) * np.cos(t) + np.sin(alpha) * np.cos(l0) * np.sin(t) + np.cos(alpha) * np.sin(b0) * np.sin(l0)
    z = np.sin(alpha) * np.sin(b0) * np.cos(t) + np.cos(alpha) * np.cos(b0)
    l = -np.arctan2(y, x)
    b = np.arcsin(z)
    return [l * 180 / np.pi, b * 180 / np.pi]


# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='moll', lon_0=0, celestial=True, resolution='c')
# draw parallels and meridians.
m.drawparallels(np.arange(-90., 90., 15.), labels=[1, 0, 0, 0], color='k')
m.drawmeridians(np.arange(0., 360., 30.))
meridian = np.arange(0., 360., 30.)
for i in np.arange(len(meridian)):
    plt.annotate(np.str(meridian[i]) + '$^\\circ$', xy=m(meridian[i], 0), xycoords='data')

l0, b0 = m(init_evs_l, init_evs_b)
m.scatter(l0, b0, s=40, c='k', marker='o', label='CRs E > 1e20 eV')

SB_l, SB_b = m(starburst[:, 0], starburst[:, 1])
m.scatter(SB_l, SB_b, s=100, marker='*', color='k', label='StarBurstGal')

AGNs_l, AGNs_b = m(AGNs[:, 0], AGNs[:, 1])
m.scatter(AGNs_l, AGNs_b, s=100, marker='^', color='m', label='AGN')

evs6_l, evs6_b = m(evs6[0], evs6[1])
m.scatter(evs6_l, evs6_b, s=20, c='b', marker='o')

from math import hypot
for i in range(len(evs1[0])):
    # __________Кола розмазування у випадковому гал полі для z=26__________
    r = evs6[2][i]
    res = circle(np.sqrt(r ** 2 + 30 ** 2), evs6[0][i], evs6[1][i]) # зсув 90 градусів вже враховано всередині функціїї
                                                                    # circ(). підставляємо просто b а не 90-b
    # Розбиваємо кожне коло на дві частини: 0 < l < 180 і l > 180
    # для того, що не малювалися лінії, які з'єднують точки
    # через всю карту
    res_left = [[], []]
    res_right = [[], []]
    for i in range(len(res[0])):
        if 0 < res[0][i] < 180:
            res_left[0].append(res[0][i])
            res_left[1].append(res[1][i])
        else:
            res_right[0].append(res[0][i])
            res_right[1].append(res[1][i])
    # Кожну з двох частин кола будемо малювати відрізками між двома сусідніми точками,
    # за умови, що відстань між ними менше певної величини. Ця величина залежить від загальної кількості
    # точок кола, яка задається всередині функції circle(): t = np.linspace(0, 2 * np.pi, 700) (коло має 700 точок)
    for i in range(len(res_left[0]) - 1):
        if hypot(res_left[0][i + 1] - res_left[0][i], res_left[1][i + 1] - res_left[1][i]) < 10:
            x_l, y_l = m(res_left[0][i], res_left[1][i])
            x_l2, y_l2 = m(res_left[0][i + 1], res_left[1][i + 1])
            m.plot([x_l, x_l2], [y_l, y_l2], color='b', ls='-', lw=1, label='z=6 + extragal')
    for i in range(len(res_right[0]) - 1):
        if hypot(res_right[0][i + 1] - res_right[0][i], res_right[1][i + 1] - res_right[1][i]) < 10:
            x_r, y_r = m(res_right[0][i], res_right[1][i])
            x_r2, y_r2 = m(res_right[0][i + 1], res_right[1][i + 1])
            m.plot([x_r, x_r2], [y_r, y_r2], color='b', ls='-', lw=1)


# _________Додавання легенди до малюнка (з видалленям дублікатів однакових ліній)_________
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
    if label not in newLabels:
        newLabels.append(label)
        newHandles.append(handle)
plt.legend(newHandles, newLabels, loc='upper right', prop={'size': 12})


plt.tight_layout()

plt.show()

