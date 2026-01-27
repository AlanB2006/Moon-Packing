def Cart2Orb(ps0,ps1,mcen):
  ps_diff = ps1 - ps0
  tempsim = rebound.Simulation()
  tempsim.units = ('AU', 'yr', 'Msun')
  tempsim.add(m=mcen)
  tempsim.add(m=0,x=ps_diff.x,y=ps_diff.y,z=ps_diff.z,vx=ps_diff.vx,vy=ps_diff.vy,vz=ps_diff.vz)
  tempsim.move_to_com()
  ps = tempsim.particles
  return ps[1].a,ps[1].e,ps[1].d,ps[1].P

earth_mass = 3.003e-6 #Solar M
r_earth = 4.259e-5
r_roche = 2.44*r_earth*(5.515/3.3)**(1/3.) #Roche Radius
a_m = 2 * r_roche
R_H = 1*(earth_mass/3)**(1/3.)
Pmoon = np.sqrt(a_m**3/earth_mass) # Kepler 3rd Law
M_moon = earth_mass/81. #Luna moon
fs = 'x-large'
tau = [1, 2, 3]

nmoon = int(input(" number of moons: "))
moon_name = input(" moon mass (L, C, P): ")

fig, axes = plt.subplots(3, 3, figsize=(12, 6), dpi=150)
for col, i in enumerate(tau):
    sa = rebound.Simulationarchive("/content/drive/MyDrive/Case-Study/Luna/Simulation Archive/Tau_698/Luna_698!_%i.bin" % i)
    npts = int(len(sa)/10)+1
    semi = np.zeros((npts,nmoon))
    ecc = np.zeros((npts,nmoon))
    p1 = np.zeros(npts)
    p2 = np.zeros(npts)
    p3 = np.zeros(npts)
    p4 = np.zeros(npts)
    p5 = np.zeros(npts)
    t = np.zeros(npts)

    cnt = 0
    for j in range(0,len(sa),10):
        ps = sa[j].particles
        if moon_name == "L":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
            ecc[cnt,0] = moon1_ecc
            semi[cnt,0] = moon1_semi/R_H
            p1[cnt] = P1

            semi[cnt,1] = moon2_semi/R_H
            ecc[cnt,1] = moon2_ecc
            p2[cnt] = P2
        if moon_name == "P":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
            moon3_semi, moon3_ecc, moon3_d, P3 = Cart2Orb(ps[1], ps[4], earth_mass + M_moon)
            moon4_semi, moon4_ecc, moon4_d, P4 = Cart2Orb(ps[1], ps[5], earth_mass + M_moon)

            ecc[cnt,0] = moon1_ecc
            semi[cnt,0] = moon1_semi/R_H
            p1[cnt] = P1

            semi[cnt,1] = moon2_semi/R_H
            ecc[cnt,1] = moon2_ecc
            p2[cnt] = P2

            semi[cnt,2] = moon3_semi/R_H
            ecc[cnt,2] = moon3_ecc
            p3[cnt] = P3

            semi[cnt,3] = moon4_semi/R_H
            ecc[cnt,3] = moon4_ecc
            p4[cnt] = P4
        if moon_name == "C":
            moon1_semi, moon1_ecc, moon1_d, P1 = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
            moon2_semi, moon2_ecc, moon2_d, P2 = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
            moon3_semi, moon3_ecc, moon3_d, P3 = Cart2Orb(ps[1], ps[4], earth_mass + M_moon)
            moon4_semi, moon4_ecc, moon4_d, P4 = Cart2Orb(ps[1], ps[5], earth_mass + M_moon)
            moon5_semi, moon5_ecc, moon5_d, P5 = Cart2Orb(ps[1], ps[6], earth_mass + M_moon)

            ecc[cnt,0] = moon1_ecc
            semi[cnt,0] = moon1_semi/R_H
            p1[cnt] = P1

            semi[cnt,1] = moon2_semi/R_H
            ecc[cnt,1] = moon2_ecc
            p2[cnt] = P2

            semi[cnt,2] = moon3_semi/R_H
            ecc[cnt,2] = moon3_ecc
            p3[cnt] = P3

            semi[cnt,3] = moon4_semi/R_H
            ecc[cnt,3] = moon4_ecc
            p4[cnt] = P4

            semi[cnt,4] = moon5_semi/R_H
            ecc[cnt,4] = moon5_ecc
            p5[cnt] = P5

        t[cnt] = sa[j].t
        cnt +=1

    color = ['r','b','g','purple','orange']
    s = .5
    ylabels = ['$a\\ (R_H)$','$e$','$p_{i+1}/p_i$']
    ylims = [(1e-2,0.45),(1e-5, 1),(1,70)]
    sublb = ['e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
    for m in range(nmoon): #moon number
        for r in range(3):
            axes[r,col].minorticks_on()
            axes[r,col].text(0.03, 0.88, sublb[r * len(tau) + col], color='k', fontsize='large', weight='bold',horizontalalignment='left', transform=axes[r,col].transAxes)
            axes[r,col].tick_params(which='major',axis='both', direction='out',length = 6.0, width = 4.0)
            axes[r,col].tick_params(which='minor',axis='both', direction='out',length = 3.0, width = 2.0)
            axes[r,col].set_ylim(ylims[r])
            axes[r,col].set_xlim(0,25)
            if r == 0:
                axes[0,col].scatter(t, semi[:, m],marker='o',color=color[m],s = s,alpha=0.5)
            elif r == 1:
                axes[1,col].scatter(t, ecc[:, m], marker='o', color=color[m], s = s, alpha=0.4)
                axes[1,col].set_yscale('log')
                locmaj = LogLocator(base=10.0, numticks=100)
                axes[1,col].yaxis.set_major_locator(locmaj)
                # 2. Define the minor ticks (2, 3, 4...9 for each decade)
                # subs=np.arange(2, 10) creates ticks at 0.2, 0.3... and 2, 3...
                locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100)
                axes[1,col].yaxis.set_minor_locator(locmin)
                axes[1,col].yaxis.set_minor_formatter(NullFormatter())
            elif r == 2:
                if m ==1:
                    axes[2, col].scatter(t, p2/p1 ,marker='o',color='b',s = s,alpha=0.5)
                elif m ==2:
                    axes[2,col].scatter(t, p3/p2 ,marker='o',color='g',s = s,alpha=0.5)
                elif m ==3:
                    axes[2,col].scatter(t, p4/p3 ,marker='o',color='purple',s = s,alpha=0.5)
                elif m ==4:
                    axes[2,col].scatter(t, p5/p4, marker='o', color='orange', s = s, alpha=0.5)
                axes[2,col].set_xlabel('Time (yr)', fontsize=fs)
            if col == 0:
                axes[r,col].set_ylabel(ylabels[r], fontsize=fs)

fig.subplots_adjust(wspace=0.3)
filename = "Luna_sa_698.png"
plt.savefig(filename, dpi=300)
plt.show()
