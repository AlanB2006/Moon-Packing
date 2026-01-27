tau = [1,2,3]
fig, axes = plt.subplots(1, 3, figsize=(12, 6), dpi=150)
earth_mass = 3.003e-6 #Solar M
R_H = 1*(earth_mass/3)**(1/3.) #0.005
color = ['r','b','g','purple','orange']
fs = 'large'
sublbl = ['b','c','d']
moon_name = input(" moon mass (L, C, P): ")

for col, i in enumerate(tau):
    ax = axes[col]
    sa = rebound.Simulationarchive("/content/drive/MyDrive/Case-Study/Pluto/Simulation Archive/Tau_698/Pluto_698_%i.bin" % i)
    m1x, m1y = [], []
    m2x, m2y = [], []
    m3x, m3y = [], []
    m4x, m4y = [], []
    m5x, m5y = [], []
    for sim in sa:
        ps = sim.particles
        if moon_name == "L":
            m1x.append((ps[2].x - ps[1].x)/R_H), m1y.append((ps[2].y - ps[1].y)/R_H)
            m2x.append((ps[3].x - ps[1].x)/R_H), m2y.append((ps[3].y - ps[1].y)/R_H)
        elif moon_name == "P":
            m1x.append((ps[2].x - ps[1].x)/R_H), m1y.append((ps[2].y - ps[1].y)/R_H)
            m2x.append((ps[3].x - ps[1].x)/R_H), m2y.append((ps[3].y - ps[1].y)/R_H)
            m3x.append((ps[4].x - ps[1].x)/R_H), m3y.append((ps[4].y - ps[1].y)/R_H)
            m4x.append((ps[5].x - ps[1].x)/R_H), m4y.append((ps[5].y - ps[1].y)/R_H)
        elif moon_name == "C":
            m1x.append((ps[2].x - ps[1].x)/R_H), m1y.append((ps[2].y - ps[1].y)/R_H)
            m2x.append((ps[3].x - ps[1].x)/R_H), m2y.append((ps[3].y - ps[1].y)/R_H)
            m3x.append((ps[4].x - ps[1].x)/R_H), m3y.append((ps[4].y - ps[1].y)/R_H)
            m4x.append((ps[5].x - ps[1].x)/R_H), m4y.append((ps[5].y - ps[1].y)/R_H)
            m5x.append((ps[6].x - ps[1].x)/R_H), m5y.append((ps[6].y - ps[1].y)/R_H)
    if moon_name == "L":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2')
        filename = "Luna_Sa.png"
    elif moon_name == "P":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2')
        ax.scatter(m3x, m3y, s=1, color='g', alpha=0.7, label='Moon 3')
        ax.scatter(m4x, m4y, s=1, color='purple', alpha=0.7, label='Moon 4')
        filename = "Pluto_Sa.png"
    elif moon_name == "C":
        ax.scatter(m1x, m1y, s=1, color='r', alpha=0.7, label='Moon 1')
        ax.scatter(m2x, m2y, s=1, color='b', alpha=0.7, label='Moon 2')
        ax.scatter(m3x, m3y, s=1, color='g', alpha=0.7, label='Moon 3')
        ax.scatter(m4x, m4y, s=1, color='purple', alpha=0.7, label='Moon 4')
        ax.scatter(m5x, m5y, s=1, color='orange', alpha=0.7, label='Moon 5')
        filename = "Ceres_Sa.png"
    else:
        print("Steelers are awesome")
    ax.set_aspect('equal')
    ax.minorticks_on()
    ax.tick_params(which='major',axis='both', direction='out',length = 6.0, width = 3.0,labelsize=fs)
    ax.tick_params(which='minor',axis='both', direction='out',length = 4.0, width = 2.0)
    ax.set_xlabel("X ($R_H$)")
    if col == 0:
        ax.set_ylabel("Y ($R_H$)")
        ax.legend(loc='lower right')
    else:
        ax.set_ylabel("")
    ax.text(0.05, 0.9, sublbl[col], color='k', fontsize=fs, weight='bold',horizontalalignment='left', transform=ax.transAxes)
    ax.set_ylim(-0.4, 0.4)
    ax.set_xlim(-0.4, 0.4)
    ax.set_aspect('equal')

plt.savefig(filename, dpi=300)
# plt.tight_layout()
plt.show()
