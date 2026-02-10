import sys, time, multiprocessing as mp, numpy as np, rebound, reboundx
from scipy.constants import astronomical_unit

lock = mp.Lock()

def write_to_file(fname, text, mode="a"):
    lock.acquire()
    f = open(fname, mode)
    f.write(text)
    f.close()
    lock.release()

# usage: python script.py N_MOONS MOON_TYPE[L/C/P] TIDES[0/1] [NPROC]
N_MOONS   = int(sys.argv[1])
MOON_TYPE = sys.argv[2]
TIDES_ON  = len(sys.argv)>3 and sys.argv[3] not in ("0","false","F","False","n","N")
TAU = float(sys.argv[4]) if (TIDES_ON and len(sys.argv) > 4) else 0.0
NPROC     = int(sys.argv[5]) if len(sys.argv)>5 else mp.cpu_count()

# embed setup in filename so multiple runs don't overwrite each other
TFLAG   = ("T1_tau%s" % TAU) if TIDES_ON else "T0"
OUTFILE = "output_%s_%dmoons_%s.txt" % (MOON_TYPE, N_MOONS, TFLAG)

e_m = 1e-7
M_S = 1.0
M_E = 3.003e-6
A_E = 1.0
e_E = 0.01671022
R_E = 6371e3/astronomical_unit
k2_E = 0.298
C_E  = 0.3308
R_H  = A_E*(M_E/3)**(1/3)

MOON_TYPES = {
    "L":{"M":M_E/81.,      "rho":3.3,  "R":1737e3/astronomical_unit},
    "C":{"M":4.72e-10,     "rho":2.2,  "R":473e3/astronomical_unit},
    "P":{"M":0.002183*M_E, "rho":1.86, "R":1.19e6/astronomical_unit}
}

MOON  = MOON_TYPES[MOON_TYPE]
M_M   = MOON["M"]
rho_M = MOON["rho"]
R_M   = MOON["R"]

R_roche = 2.44*R_E*(5.515/rho_M)**(1/3)
A_M0    = 2*R_roche
# reference period at 1.8 * R_roche (years), used for timestep and check spacing
P_REF = np.sqrt((1.8*R_roche)**3 / M_E)

def cart2orb(ps0, ps1, mcen):
    ps_diff = ps1-ps0
    s = rebound.Simulation()
    s.units=("AU","yr","Msun")
    s.add(m=mcen)
    s.add(m=0,x=ps_diff.x,y=ps_diff.y,z=ps_diff.z,vx=ps_diff.vx,vy=ps_diff.vy,vz=ps_diff.vz)
    s.move_to_com()
    p = s.particles
    return p[1].a, p[1].e, p[1].d

def get_semi(dH,a_j,k,j):
    X=0.5*((2*M_M)/(3*(M_E+k*M_M)))**(1/3)
    return a_j*((1+dH*X)/(1-dH*X))**(k-j)

def get_true_anomaly(j):
    phi=(1+np.sqrt(5))/2
    return (j*(phi*2*np.pi))%(2*np.pi)

def setup_tides(sim):
    ps=sim.particles
    rbx=reboundx.Extras(sim)
    ts=rbx.load_force("tides_spin")
    rbx.add_force(ts)
    ps["Earth"].params["k2"]=k2_E/2
    ps["Earth"].params["I"]=C_E*M_E*R_E**2
    spin_period=6/24/365.25
    spin_mag=2*np.pi/spin_period
    obl=np.radians(23.44)
    ps["Earth"].params["Omega"]=[0, spin_mag*np.sin(obl), spin_mag*np.cos(obl)]
    ps["Earth"].params["tau"]=TAU/3600/24/365.25
    L=np.array(sim.angular_momentum())+np.array(rbx.spin_angular_momentum())
    rot=rebound.Rotation.to_new_axes(newz=L)
    rbx.rotate_simulation(rot)
    rbx.initialize_spin_ode(ts)

def run_simulation(delta_H):
    t0=time.time()
    # delta_H = delta_rng[d_idx]
    sim = rebound.Simulation()
    sim.collision = "direct"
    sim.collision_resolve = "halt"
    sim.units = ('AU', 'yr', 'Msun')
    sim.integrator = "ias15"
    sim.ri_trace.peri_mode = 1

    am  = np.array([A_M0]+[get_semi(delta_H,A_M0,i,0) for i in range(1,N_MOONS)])
    am0 = am.copy()  # store initial semimajor axes

    sim.add(m=M_S,hash="Sun")
    sim.add(m=M_E,a=A_E,e=e_E,r=R_E,primary=sim.particles[0],hash="Earth")
    for i in range(N_MOONS):
        sim.add(m=M_M,a=am[i],e=e_m,r=R_M,f=get_true_anomaly(i),primary=sim.particles["Earth"],hash="M%d"%(i+1))
    sim.move_to_com()

    if TIDES_ON: setup_tides(sim)
    ps=sim.particles

    P_M=np.sqrt(am[0]**3/M_E)          # period at actual inner moon (2*R_roche)
    tmax=1e7*P_M                       # total run length in inner-moon orbits
    sim.dt=P_REF/20.0                  # timestep from 1.8*R_roche period
    sim.ri_ias15.min_dt=P_REF/50.0
    dt_far,dt_near=1000*P_REF,40*P_REF # check spacing based on P_REF
    dt_check=dt_far

    emax=np.zeros(N_MOONS)
    emin=np.ones(N_MOONS)
    t=0.0 # Start time
    fate_label="stable"

    try:
        while t<tmax:
            t_target=t+dt_check
            if t_target>tmax: t_target=tmax
            sim.integrate(t_target)
            t=sim.t
            close=False
            for i in range(N_MOONS):
                a_i,e_i,d_i=cart2orb(ps["Earth"],ps[2+i],M_E+M_M)
                if e_i>emax[i]: emax[i]=e_i
                if e_i<emin[i]: emin[i]=e_i
                if a_i<0:       fate_label="a<0";       t=tmax; break
                if e_i>1:       fate_label="e>1";       t=tmax; break
                if d_i<R_roche: fate_label="d<R_roche"; t=tmax; break
                if a_i>0.4*R_H: fate_label="a>0.4R_H";  t=tmax; break
                if d_i<3*R_roche or a_i>0.35*R_H: close=True
            if t>=tmax: break
            dt_check = dt_near if close else dt_far
        if fate_label=="stable" and t>=tmax: fate_label="tmax"
    except rebound.Collision:
        fate_label="collision"

    wall=time.time()-t0
    out = "%1.2f, " % (delta_H,)
    for i in range(N_MOONS):
        out += "%1.4f, " % (emax[i])
    out += "%1.4e, %1.5f, %s\n" % (sim.t/P_M, wall, fate_label)
    write_to_file(OUTFILE, out)
    return 1

# header
hdr = "# dH, "
for i in range(N_MOONS):
    hdr += "eM%d_max, " % (i+1)
hdr += "t_end(PM), wall_time, fate\n"
write_to_file(OUTFILE, hdr, mode="w")

# delta_H
par_dH = np.arange(3.5,14,0.01) #Luna: 9, Pluto: 7.5, Ceres: 14.0
PARAMS = par_dH
N_TOTAL=len(PARAMS)

print("%d runs (%d moons, %s, tides=%s) using %d procs → %s"
      % (N_TOTAL, N_MOONS, MOON_TYPE, str(TIDES_ON), NPROC, OUTFILE))
with mp.Pool(processes=NPROC) as pool:
    done=0; BAR=30
    for _ in pool.imap_unordered(run_simulation, PARAMS):
        done+=1
        frac=done/float(N_TOTAL)
        bar="#"*int(frac*BAR)+"-"*(BAR-int(frac*BAR))
        print("\r[%s] %d/%d (%5.1f%%)" % (bar, done, N_TOTAL, 100*frac), end="")
print("\nDone.")
