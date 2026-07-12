import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt
import threading
import os
from scipy.constants import astronomical_unit

#https://github.com/dtamayo/reboundx/blob/main/ipython_examples/TidesSpinEarthMoon.ipynb

def write2file(fname,f_str,head):
    lock.acquire() # thread blocks at this line until it can obtain lock
    if head == "head":
            f = open(fname,'w')
    else:
        f = open(fname,'a')
    f.write(f_str)
    lock.release()

def Cart2Orb(ps0,ps1,Mcen):
    temp_sim = rebound.Simulation()
    temp_sim.units = ('AU','Msun','years')
    temp_sim.add(m=Mcen,x=ps0.x,y=ps0.y,z=ps0.z,vx=ps0.vx,vy=ps0.vy,vz=ps0.vz)
    temp_sim.add(x=ps1.x,y=ps1.y,z=ps1.z,vx=ps1.vx,vy=ps1.vy,vz=ps1.vz)
    temp_sim.move_to_hel()
    temp_ps = temp_sim.particles

    return temp_ps[1].a, temp_ps[1].P*365.25, temp_ps[1].e

fname = "reb_tide.txt"
lock = threading.Lock()
write2file(fname,"",'head')

R_sun = 0.00465047 #radius of Sun in AU

# Alpha Cen parameters; https://ui.adsabs.harvard.edu/abs/2021AJ....162...14A/abstract
M_A = 1.0 # M_sun
R_A = R_sun

#tidal + earth parameters
m_earth = 3e-6
r_earth = 4.259e-5
k2_earth = 0.298
c_earth = 0.3308 # Dimensionless moment of inertia, normalized by MR^2
a_earth = 1
P_earth = np.sqrt(a_earth**3/M_A)
spin_period_earth = (6./ 24. / 365.25) # 1 day in default REBOUND units
spin_mag_earth = 2 * np.pi / spin_period_earth
obl_earth = np.radians(23.44)

#tidal + moon parameters
r_roche = 2.44*r_earth*(5.515/3.3)**(1/3.) #Roche Radius 1.233e-04 #Luna 3.3 Ceres 2.2 Pluto 1.86
m_moon = m_earth/81.
r_moon = 1737e3/astronomical_unit
k2_moon = 0.239
c_moon = 0.394
a_moon = 2 * r_roche  #5*r_earth
P_moon = np.sqrt(a_moon**3/(m_earth+m_moon))
MA_moon = np.pi/3.
R_H = 1*(m_earth/3)**(1/3.)


sim = rebound.Simulation()
sim.units = ('AU','Msun','years')
sim.integrator = 'ias15'
#sim.ri_trace.peri_mode = "FULL_IAS15"
sim.dt = P_moon/5.

sim.add(m=M_A,hash='Star')
sim.add(m=m_earth,r=r_earth,a=a_earth,e=0.0167,primary=sim.particles[0],hash='Earth')
sim.add(m=m_moon,a=a_moon,e=0.0001,M=MA_moon,primary=sim.particles[1],hash='Moon')

sim.move_to_com()
ps = sim.particles

rebx = reboundx.Extras(sim)
ts = rebx.load_force("tides_spin")
rebx.add_force(ts)

ps["Earth"].params['k2'] = k2_earth/2
ps["Earth"].params['I'] = c_earth * m_earth * r_earth**2 # if m=0, the effect interprets 'I' as I/m, so we could set ...['I'] = c_earth * r_earth**2
ps["Earth"].params['Omega'] = [0., spin_mag_earth * np.sin(obl_earth), spin_mag_earth * np.cos(obl_earth)]
ps["Earth"].params['tau'] = 698./3600./24./365.25

Ltot = np.array(sim.angular_momentum()) + np.array(rebx.spin_angular_momentum())
rot = rebound.Rotation.to_new_axes(newz=Ltot) # Alignment to invariant plane
rebx.rotate_simulation(rot)
rebx.initialize_spin_ode(ts)

t_fin = 1e5 # REBOUND time units
#times = np.concatenate((np.arange(0,1000,100),np.arange(1000,t_fin+1000,1000)))
times = np.arange(0,3e4,100)
Nout = len(times)
lunar_a = np.zeros(Nout)
lunar_orbit_period = np.zeros(Nout)

earth_sv_hat = np.zeros((Nout, 3))
earth_spin_mag = np.zeros(Nout)
earth_obliquity = np.zeros(Nout)
earth_phase = np.zeros(Nout)

earth_sv = np.zeros((Nout, 3))

for i, t in enumerate(times):
    earth_spin_axis_inv = ps["Earth"].params['Omega'] # Initially invarant frame
    earth_orbit_normal = ps["Earth"].hvec
    rot_to_orbit_frame = rebound.Rotation.to_new_axes(newz=earth_orbit_normal)
    earth_spin_axis_orb = rot_to_orbit_frame * earth_spin_axis_inv # Now in the Earth's orbit frame

    # We save the unit earth spin vector
    earth_sv_hat[i] = earth_spin_axis_orb / np.linalg.norm(earth_spin_axis_orb)
    # We can also interpret the spin axis in the more natural spherical coordinates
    mag, obliquity, phase = rebound.xyz_to_spherical(earth_spin_axis_orb)

    earth_spin_mag[i] = mag
    earth_obliquity[i] = np.degrees(obliquity)
    earth_phase[i] = np.degrees(phase)

    sim.integrate(t)
    # Orbital Parameters
    lunar_a[i], lunar_orbit_period[i], e_m = Cart2Orb(ps[1],ps[2],m_earth+m_moon)
    out_stg = "%1.5e,%1.4f,%1.4f,%1.4f,%1.4f,%1.6e\n" % (t,lunar_a[i]/R_H,e_m,ps[1].a,ps[1].e,earth_spin_mag[i])
    write2file(fname,out_stg,'foot')