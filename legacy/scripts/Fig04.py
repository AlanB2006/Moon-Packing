import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogLocator, NullFormatter
import rebound

def Cart2Orb(ps0, ps1, mcen):
    ps_diff = ps1 - ps0
    tempsim = rebound.Simulation()
    tempsim.units = ("AU", "yr", "Msun")
    tempsim.add(m=mcen)
    tempsim.add(
        m=0,
        x=ps_diff.x, y=ps_diff.y, z=ps_diff.z,
        vx=ps_diff.vx, vy=ps_diff.vy, vz=ps_diff.vz
    )
    tempsim.move_to_com()
    ps = tempsim.particles
    return ps[1].a, ps[1].e, ps[1].d, ps[1].P

tau = [1, 2, 3]
nmoon = int(input(" number of moons: "))
moon_name = input(" moon mass (L, C, P): ").upper()

earth_mass = 3.003e-6  # Solar M
r_earth = 4.259e-5
r_roche = 2.44 * r_earth * (5.515 / 3.3) ** (1 / 3.0)
a_m = 2 * r_roche
R_H = 1 * (earth_mass / 3) ** (1 / 3.0)
Pmoon = np.sqrt(a_m**3 / earth_mass)
M_moon = earth_mass / 81.0

fig = plt.figure(figsize=(12, 13), dpi=150)

gs = GridSpec(
    4, 3, figure=fig,
    height_ratios=[2.2, 1, 1, 1],
    hspace=0.15, wspace=0.15
)

axes_top = [fig.add_subplot(gs[0, c]) for c in range(3)]
axes_bottom = np.array(
    [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(1, 4)]
)
s = 0.5
fs_bottom = 'x-large'
color = ['r', 'b', 'brownp', 'purple', 'orange']
fs_top = 'large'
sublbl_top = ['a', 'b', 'c']

axes = axes_top

for col, i in enumerate(tau):
    ax = axes[col]

    sa = rebound.Simulationarchive(
        "/content/C_%i.bin" % i
    )

    m1x, m1y = [], []
    m2x, m2y = [], []
    m3x, m3y = [], []
    m4x, m4y = [], []
    m5x, m5y = [], []

    for sim in sa:
        ps = sim.particles
        if moon_name == "L":
            m1x.append((ps[2].x - ps[1].x) / R_H); m1y.append((ps[2].y - ps[1].y) / R_H)
            m2x.append((ps[3].x - ps[1].x) / R_H); m2y.append((ps[3].y - ps[1].y) / R_H)
        elif moon_name == "P":
            m1x.append((ps[2].x - ps[1].x) / R_H); m1y.append((ps[2].y - ps[1].y) / R_H)
            m2x.append((ps[3].x - ps[1].x) / R_H); m2y.append((ps[3].y - ps[1].y) / R_H)
            m3x.append((ps[4].x - ps[1].x) / R_H); m3y.append((ps[4].y - ps[1].y) / R_H)
        elif moon_name == "C":
            m1x.append((ps[2].x - ps[1].x) / R_H); m1y.append((ps[2].y - ps[1].y) / R_H)
            m2x.append((ps[3].x - ps[1].x) / R_H); m2y.append((ps[3].y - ps[1].y) / R_H)
            m3x.append((ps[4].x - ps[1].x) / R_H); m3y.append((ps[4].y - ps[1].y) / R_H)
            m4x.append((ps[5].x - ps[1].x) / R_H); m4y.append((ps[5].y - ps[1].y) / R_H)
            m5x.append((ps[6].x - ps[1].x) / R_H); m5y.append((ps[6].y - ps[1].y) / R_H)

    # axins = ax.inset_axes([0.15, 0.1, 0.3, 0.3]) #left bottom width heght

    if moon_name == "L":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2') #main fig
        # axins.scatter(m1x, m1y, s=1, color='r', alpha=0.7) #zoom fig
        # axins.scatter(m2x, m2y, s=1, color='b', alpha=0.7)
    elif moon_name == "P":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2')
        ax.scatter(m3x, m3y, s=1, color='saddlebrown', alpha=0.7, label='Moon 3')
        # axins.scatter(m1x, m1y, s=1, color='r', alpha=0.7)
        # axins.scatter(m2x, m2y, s=1, color='b', alpha=0.7)
        # axins.scatter(m3x, m3y, s=1, color='saddlebrown', alpha=0.7)
    elif moon_name == "C":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2')
        ax.scatter(m3x, m3y, s=1, color='saddlebrown', alpha=0.7, label='Moon 3')
        ax.scatter(m4x, m4y, s=1, color='purple', alpha=0.7, label='Moon 4')
        ax.scatter(m5x, m5y, s=1, color='orange', alpha=0.7, label='Moon 5')
        # axins.scatter(m1x, m1y, s=1, color='r', alpha=0.7)
        # axins.scatter(m2x, m2y, s=1, color='b', alpha=0.7)
        # axins.scatter(m3x, m3y, s=1, color='saddlebrown', alpha=0.7)
        # axins.scatter(m4x, m4y, s=1, color='purple', alpha=0.7)
        # axins.scatter(m5x, m5y, s=1, color='orange', alpha=0.7)
    else:
        raise ValueError("moon_name must be one of: L, P, C")

    # ax.set_aspect('auto')
    ax.minorticks_on()
    ax.tick_params(which='major', axis='both', direction='out', length=6.0, width=4.0)
    ax.tick_params(which='minor', axis='both', direction='out', length=3.0, width=2.0)
    ax.set_xlabel("X ($R_H$)", fontsize=fs_bottom)


    if col != 0:
        ax.tick_params(labelleft=False)
    if col == 0:
        ax.set_ylabel("Y ($R_H$)", fontsize=fs_bottom)
        ax.legend(loc='upper right', markerscale=5)
        # axins = ax.inset_axes([0.67, 0.085, 0.3, 0.3]) #left bottom width heght
        # axins.scatter(m1x, m1y, s=1, color='r', alpha=0.7)
        # axins.scatter(m2x, m2y, s=1, color='b', alpha=0.7)
        # axins.scatter(m3x, m3y, s=1, color='saddlebrown', alpha=0.7)
        # axins.scatter(m4x, m4y, s=1, color='purple', alpha=0.7)
        # axins.scatter(m5x, m5y, s=1, color='orange', alpha=0.7)
        # axins.set_xlim(-0.1, 0.1)
        # axins.set_ylim(-0.1, 0.1)
        # axins.tick_params(labelsize=8)
    else:
        ax.set_ylabel("")

    ax.text(0.05, 0.9, sublbl_top[col], color='k', fontsize=fs_top, weight='bold',
            horizontalalignment='left', transform=ax.transAxes)
    ax.set_ylim(-0.44, 0.44)
    ax.set_xlim(-0.44, 0.44)

    ax.set_box_aspect(1)
    ax.set_aspect('equal', adjustable='box')

color = ['r', 'b', 'saddlebrown', 'purple', 'orange']
ylabels = ['$a\ (R_H)$', '$e$', '$\\frac{p_{i+1}}{p_i}$']
ylims = [(0, 0.42), (1e-3, 1), (1.3, 2.9)]
sublb = ['d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

axes = axes_bottom

for col, i in enumerate(tau):
    sa = rebound.Simulationarchive(
        "/content/C_%i.bin" % i)

    npts = int(len(sa) / 10) + 1
    semi = np.zeros((npts, nmoon))
    ecc = np.zeros((npts, nmoon))
    p1 = np.zeros(npts)
    p2 = np.zeros(npts)
    p3 = np.zeros(npts)
    p4 = np.zeros(npts)
    p5 = np.zeros(npts)
    t = np.zeros(npts)

    cnt = 0
    for j in range(0, len(sa), 10):
        ps = sa[j].particles

        if moon_name == "L":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)

            ecc[cnt, 0] = moon1_ecc
            semi[cnt, 0] = moon1_semi / R_H
            p1[cnt] = P1

            semi[cnt, 1] = moon2_semi / R_H
            ecc[cnt, 1] = moon2_ecc
            p2[cnt] = P2

        if moon_name == "P":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
            moon3_semi, moon3_ecc, moon3_d, P3 = Cart2Orb(ps[1], ps[4], earth_mass + M_moon)

            ecc[cnt, 0] = moon1_ecc
            semi[cnt, 0] = moon1_semi / R_H
            p1[cnt] = P1

            semi[cnt, 1] = moon2_semi / R_H
            ecc[cnt, 1] = moon2_ecc
            p2[cnt] = P2

            semi[cnt, 2] = moon3_semi / R_H
            ecc[cnt, 2] = moon3_ecc
            p3[cnt] = P3

        if moon_name == "C":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
            moon3_semi, moon3_ecc, moon3_d, P3 = Cart2Orb(ps[1], ps[4], earth_mass + M_moon)
            moon4_semi, moon4_ecc, moon4_d, P4 = Cart2Orb(ps[1], ps[5], earth_mass + M_moon)
            moon5_semi, moon5_ecc, moon5_d, P5 = Cart2Orb(ps[1], ps[6], earth_mass + M_moon)

            ecc[cnt, 0] = moon1_ecc
            semi[cnt, 0] = moon1_semi / R_H
            p1[cnt] = P1

            semi[cnt, 1] = moon2_semi / R_H
            ecc[cnt, 1] = moon2_ecc
            p2[cnt] = P2

            semi[cnt, 2] = moon3_semi / R_H
            ecc[cnt, 2] = moon3_ecc
            p3[cnt] = P3

            semi[cnt, 3] = moon4_semi / R_H
            ecc[cnt, 3] = moon4_ecc
            p4[cnt] = P4

            semi[cnt, 4] = moon5_semi / R_H
            ecc[cnt, 4] = moon5_ecc
            p5[cnt] = P5

        t[cnt] = sa[j].t
        cnt += 1

    # plotting
    for m in range(nmoon):  # moon number
        for r in range(3):
            axes[r, col].minorticks_on()
            current_label = sublb[r * len(tau) + col]
            if current_label in [ 'l']:
                x_pos, y_pos = 0.03, 0.88#0.91, 0.89
                align = 'left'
            elif current_label in ["f"]:
                x_pos, y_pos = 0.03, 0.88#0.91, 0.91
                align = 'left'
            else:
                x_pos, y_pos = 0.03, 0.88
                align = 'left'
            axes[r, col].text(
                x_pos, y_pos, current_label,
                color='k', fontsize='large', weight='bold',
                horizontalalignment=align,
                transform=axes[r, col].transAxes
            )
            axes[r, col].tick_params(which='major', axis='both', direction='out', length=6.0, width=4.0)
            axes[r, col].tick_params(which='minor', axis='both', direction='out', length=3.0, width=2.0)
            axes[r, col].set_ylim(ylims[r])
            axes[r, col].set_xlim(0, 25)
            if col != 0:
                axes[r,col].tick_params(labelleft=False)
            if r != 2:
                axes[r,col].tick_params(labelbottom=False)
            if r == 0:
                axes[0, col].scatter(t, semi[:, m], marker='o', color=color[m], s=s, alpha=0.5)
                # axes[0, col].set_yscale('log')
            # if col == 0 and r == 0:
            #     axesins = axes[r, col].inset_axes([0.45, 0.45, 0.50, 0.50])

            #     for mm in range(nmoon):
            #         axesins.scatter(t, semi[:, mm], marker='o', color=color[mm], s=s, alpha=0.5)

            #     axesins.set_xlim(0, 10)          # zoomed time range
            #     axesins.set_ylim(0.0, 0.1)    # zoomed semimajoraxis range
            #     axesins.tick_params(labelsize=7)

            elif r == 1:
                axes[1, col].scatter(t, ecc[:, m], marker='o', color=color[m], s=s, alpha=0.7 if m == 2 else 0.4)
                axes[1, col].set_yscale('log')
                locmaj = LogLocator(base=10.0, numticks=100)
                axes[1, col].yaxis.set_major_locator(locmaj)
                locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100)
                axes[1, col].yaxis.set_minor_locator(locmin)
                axes[1, col].yaxis.set_minor_formatter(NullFormatter())
            elif r == 2:
                if m == 1:
                    axes[2, col].scatter(t, p2 / p1, marker='o', color='b', s=s, alpha=0.6) # make k when Luna case make b on others
                elif m == 2:
                    axes[2, col].scatter(t, p3 / p2, marker='o', color='saddlebrown', s=s, alpha=0.5)
                elif m == 3:
                    axes[2, col].scatter(t, p4 / p3, marker='o', color='purple', s=s, alpha=0.3)
                elif m == 4:
                    axes[2, col].scatter(t, p5 / p4, marker='o', color='orange', s=s, alpha=0.5)

                axes[2, col].set_xlabel('Time (yr)', fontsize=fs_bottom)

            if col == 0:
                axes[r, col].set_ylabel(ylabels[r], fontsize=fs_bottom)

            # if col == 0 and r == 2:
            #     axesins = axes[r, col].inset_axes([0.45, 0.45, 0.50, 0.50])
            #     if nmoon >= 2:
            #         axesins.scatter(t, p2/p1, marker='o', color='b', s=s, alpha=0.6)
            #     if nmoon >= 3:
            #         axesins.scatter(t, p3/p2, marker='o', color='saddlebrown', s=s, alpha=0.5)
            #     if nmoon >= 4:
            #         axesins.scatter(t, p4/p3, marker='o', color='purple', s=s, alpha=0.4)
            #     if nmoon >= 5:
            #         axesins.scatter(t, p5/p4, marker='o', color='orange', s=s, alpha=0.5)

                # axesins.set_xlim(0, 10)          # zoomed time range
                # axesins.set_ylim(1.5, 1.6)    # zoomed semimajor-axis range
                # axesins.tick_params(labelsize=7)
fig.subplots_adjust(wspace=0.05, hspace=0.25, right=0.90)
fig.savefig("Ceres_Time_Series & (x,y)", dpi=300, bbox_inches="tight")
plt.show()