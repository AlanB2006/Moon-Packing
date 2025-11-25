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
fig, axes = plt.subplots(3, 3, figsize=(12, 6), dpi=150)
for col, i in enumerate(tau):
    sa = rebound.Simulationarchive("" % i)
    npts = int(len(sa)/10)+1
    semi = np.zeros((npts,2))
    ecc = np.zeros((npts,2))
    p1 = np.zeros(npts)
    p2 = np.zeros(npts)
    t = np.zeros(npts)

    cnt = 0
    for j in range(0,len(sa),10):
        ps = sa[j].particles
        moon1_semi, moon1_ecc, moon1_d, moon_p = Cart2Orb(ps[1], ps[2], earth_mass + M_moon)
        moon2_semi, moon2_ecc, moon2_d, moon_P = Cart2Orb(ps[1], ps[3], earth_mass + M_moon)
        ecc[cnt,0] = moon1_ecc
        semi[cnt,0] = moon1_semi/R_H
        p1[cnt] = moon_p
        semi[cnt,1] = moon2_semi/R_H
        ecc[cnt,1] = moon2_ecc
        p2[cnt] = moon_P
        t[cnt] = sa[j].t
        cnt +=1

    color = ['r','b']
    ylabels = ['$a\\ (R_H)$','$e$','$p_2/p_1$']
    ylims = [(1e-2,0.45),(1e-5, 1),(1,70)]
    for m in range(2): #moon number
        for r in range(3):
            axes[r,col].minorticks_on()
            axes[r,col].tick_params(which='major',axis='both', direction='out',length = 6.0, width = 4.0)
            axes[r,col].tick_params(which='minor',axis='both', direction='out',length = 3.0, width = 2.0)
            axes[r,col].set_ylim(ylims[r])
            if r == 0:
                axes[0,col].scatter(t, semi[:, m],marker='o',color=color[m],s=1,alpha=0.5)
            elif r == 1:
                axes[1,col].scatter(t, ecc[:, m], marker='o', color=color[m], s=1, alpha=0.5)
                axes[1,col].set_yscale('log')
                # axes[1,col].yaxis.set_minor_locator(tck.LogLocator(subs=[1,2], numticks=99))
                axes[1, col].yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1.0, 2.0)*0.1, numticks=10))
                axes[1,col].yaxis.set_minor_formatter(tck.NullFormatter())
                # axes[1,col].tick_params(axis='y', which='minor', length=4, color='k', width = 0.5)
                axes[1,col].set_xlabel('Time', fontsize=fs)
            elif r == 2:
                axes[2, col].scatter(t, p2/p1,marker='o',color='black',s=1,alpha=0.5)     
            if col == 0:
                axes[r,col].set_ylabel(ylabels[r], fontsize=fs)


fig.subplots_adjust(wspace=0.3)
plt.show()
# moon_name = input(" moon mass (L, C, P): ")
# if moon_name == "L":
#     filename = "Luna_Sa.png"
# elif moon_name == "C":
#     filename = "Ceres_Sa.png"
# elif moon_name == "P":
#     filename = "Pluto_Sa.png"
# else:
#     exit()
# plt.savefig(filename, dpi=300)
