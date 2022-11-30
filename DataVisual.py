from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Events import events_1020_no_calib, events_1020_ta095_pao087, events_1020_ta1052_pao0948


def read_events(filename):
    file = open(filename, 'r')
    l = []
    b = []
    theta = []
    for line in file.readlines():
        l.append(float(line.split()[1]))
        b.append(float(line.split()[2]))
        theta.append(float(line.split()[3]))
    return np.array([l, b, theta])


init_evs_l = np.array([events_1020_ta095_pao087[i][1] for i in range(len(events_1020_ta095_pao087))])
init_evs_b = np.array([events_1020_ta095_pao087[i][2] for i in range(len(events_1020_ta095_pao087))])

evs1 = read_events('z1_events_1020_ta095_pao087_theta.txt')
evs6 = read_events('z6_events_1020_ta095_pao087_theta.txt')
evs26 = read_events('z26_events_1020_ta095_pao087_theta.txt')


''' RA, Dec are arrays of the same length.
RA takes values in [0,360), Dec in [-90,90],
which represent angles in degrees.
org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
title is the title of the figure.
projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
'''
org = 180
projection = 'mollweide'

x0 = np.remainder(init_evs_l+360-org, 360) # shift RA values
ind0 = x0 > 180
x0[ind0] -= 360    # scale conversion to [-180, 180]
x0 = -x0    # reverse the scale: East to the left

x1 = np.remainder(evs1[0]+360-org, 360) # shift RA values
ind1 = x1 > 180
x1[ind1] -= 360    # scale conversion to [-180, 180]
x1 = -x1    # reverse the scale: East to the left

x6 = np.remainder(evs6[0]+360-org, 360) # shift RA values
ind6 = x6 > 180
x6[ind6] -= 360    # scale conversion to [-180, 180]
x6 = -x6    # reverse the scale: East to the left

x26 = np.remainder(evs26[0]+360-org, 360) # shift RA values
ind26 = x26 > 180
x26[ind26] -= 360    # scale conversion to [-180, 180]
x26 = -x26    # reverse the scale: East to the left

tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210]) # signature of longitude
tick_labels = np.remainder(tick_labels+360+org, 360)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection=projection)

ax.scatter(np.radians(x0), np.radians(init_evs_b), s=30, c='k', marker='o')
ax.scatter(np.radians(x1), np.radians(evs1[1]), s=30, c='k', marker='x')
ax.scatter(np.radians(x6), np.radians(evs6[1]), s=350, c='k', marker='+')
ax.scatter(np.radians(x26), np.radians(evs26[1]), s=30, c='k', marker='D')

theta_circ = np.linspace(0, 2*np.pi, 100)

for i in range(len(evs1[0])):
    ax.annotate(s='', xy=(np.radians(x0)[i], np.radians(init_evs_b)[i]),
                xytext=(np.radians(x1)[i], np.radians(evs1[1])[i]),
                arrowprops=dict(arrowstyle='<-', shrinkA=0, shrinkB=0, linewidth=2, alpha=0.6))
    ax.annotate(s='', xy=(np.radians(x1)[i], np.radians(evs1[1])[i]),
                xytext=(np.radians(x6)[i], np.radians(evs6[1])[i]),
                arrowprops=dict(arrowstyle='<-', shrinkA=0, shrinkB=0, linewidth=2, alpha=0.6))
    ax.annotate(s='', xy=(np.radians(x6)[i], np.radians(evs6[1])[i]),
                xytext=(np.radians(x26)[i], np.radians(evs26[1])[i]),
                arrowprops=dict(arrowstyle='<-', shrinkA=0, shrinkB=0, linewidth=2, alpha=0.6))
    r = np.radians(evs26[2])[i]
    x = r * np.cos(theta_circ) - np.radians(x26)[i]
    y = r * np.sin(theta_circ) + np.radians(evs26[1])[i]
    ax.plot(-x, y, color='k', linestyle=':')
    r = np.radians(evs6[2])[i]
    x = r * np.cos(theta_circ) - np.radians(x6)[i]
    y = r * np.sin(theta_circ) + np.radians(evs6[1])[i]
    ax.plot(-x, y, color='k', linestyle='--')
    r = np.radians(evs1[2])[i]
    x = r * np.cos(theta_circ) - np.radians(x1)[i]
    y = r * np.sin(theta_circ) + np.radians(evs1[1])[i]
    ax.plot(-x, y, color='k')


ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
ax.set_title('z1-6-26_events_1020_ta095_pao087')
ax.title.set_fontsize(15)
ax.xaxis.label.set_fontsize(12)
ax.yaxis.label.set_fontsize(12)
ax.grid(color='tab:gray', linestyle='-', linewidth=0.2)
plt.show()



#plot_mwd(evs1[0], evs1[1], org=0, title='z1_events_1020_no_calib', projection='mollweide')



