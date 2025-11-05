n_moons = int(input("number of moons: "))
moon_mass = input(" moon mass (L, C, P): ")

def write_to_file(fname,f_str,head):
    lock.acquire()
    if head == "head":
      f = open(fname,"w")
    else:
      f = open(fname,"a")
    f.write(f_str)
    lock.release()

def Cart2Orb(ps0,ps1,mcen):
  ps_diff = ps1 - ps0
  tempsim = rebound.Simulation()
  tempsim.units = ('AU', 'yr', 'Msun')
  tempsim.add(m=mcen)
  tempsim.add(m=0,x=ps_diff.x,y=ps_diff.y,z=ps_diff.z,vx=ps_diff.vx,vy=ps_diff.vy,vz=ps_diff.vz)
  tempsim.move_to_com()
  ps = tempsim.particles
  return ps[1].a,ps[1].e,ps[1].d

def get_semi(delta_H,a_j,k,j):
  X = 0.5*((2*M_moon)/(3*(earth_mass+k*M_moon)))**(1./3.)
  return a_j*((1+delta_H*X)/(1-delta_H*X))**(k-j)

def get_TrueAnomaly(j):
  phi = (1+np.sqrt(5))/2.
  golden_angle =  phi * 2 * np.pi
  return (j * golden_angle) % (2 * np.pi)

collision_happened = False
def collision_callback(sim_pointer, collision):
    global collision_happened
    collision_happened = True
    return 0

def calc_betamax(N,a_N,a_1,m_p,m_i):
  mass_fact = (12*m_p/m_i)**(1/3)
  top = (a_N/a_1)**(1/(N-1)) - 1
  bot = (a_N/a_1)**(1/(N-1)) + 1
  return mass_fact*(top/bot)

def run_simulation(d_idx):
  global collision_happened
  collision_happened = False
  start_wall_time = time.time()
  delta_H = delta_rng[d_idx]
  sim = rebound.Simulation()
  sim.collision = "direct"
  sim.collision_resolve = "halt"
  sim.units = ('AU', 'yr', 'Msun')
  sim.integrator = "ias15"
  sim.ri_trace.peri_mode = 1
  am = np.array([a_m] + [get_semi(delta_H, a_m, i, 0) for i in range(1, n_moons)])
  sim.add(m=star_mass,hash='Sun') # 0
  sim.add(m=earth_mass, a=earth_semi, e=earth_ecc,r=r_earth,primary=sim.particles[0],hash='Earth') # 1
  for i in range(0,n_moons):
    sim.add(m=M_moon, a=am[i], e=e_m,r=radius_moon, f=get_TrueAnomaly(i),primary=sim.particles['Earth'],hash='Moon %i' % (i+1)) # 2 3 4
  sim.move_to_com()

  rebx = reboundx.Extras(sim)
  ts = rebx.load_force("tides_spin")
  rebx.add_force(ts)
  ps = sim.particles
  ps["Earth"].params['k2'] = k2_earth
  ps["Earth"].params['I'] = c_earth * earth_mass * r_earth**2 # if m=0, the effect interprets 'I' as I/m, so we could set ...['I'] = c_earth * r_earth**2
  spin_period_earth = (6/24. / 365.25)  # 6 hour spin period
  spin_mag_earth = 2 * np.pi / spin_period_earth
  obl_earth = np.radians(23.44)
  ps["Earth"].params['Omega'] = [0., spin_mag_earth * np.sin(obl_earth), spin_mag_earth * np.cos(obl_earth)]
  ps["Earth"].params['tau'] = 100./3600./24./365.25
  Ltot = np.array(sim.angular_momentum()) + np.array(rebx.spin_angular_momentum())
  rot = rebound.Rotation.to_new_axes(newz=Ltot) # Alignment to invariant plane
  rebx.rotate_simulation(rot)
  rebx.initialize_spin_ode(ts)

  Pmoon = np.sqrt(am[0]**3/earth_mass) # Kepler 3rd Law
  max_time = 1e7*Pmoon # Total spins
  sim.dt = Pmoon / 20  # Timestep calc per orbit
  sim.ri_ias15.min_dt = Pmoon / 50
  stable = True
  out_step = 1000*Pmoon # check every nth orbits
  emax = np.zeros(n_moons)
  t = sim.t
  sim.save_to_file("archive.bin", interval=100.)
  try:
    while stable:
      sim.integrate(t)
      current_wall_time = (time.time() - start_wall_time)
      if t >= max_time:
          stable = False
      for i in range(0,n_moons):
        moon_semi, moon_ecc, moon_d = Cart2Orb(ps['Earth'],ps[2+i],earth_mass+M_moon) # returns/inputs 
        fname = "evolution%i.txt" % i
        with open(fname, "a") as outfile:
          outstg = "%1.2f, %1.4f, %1.4f, %1.4f\n" % (sim.t, moon_semi/R_H, moon_ecc, moon_d/R_H)
          outfile.write(outstg)
        if moon_ecc > emax[i]:  emax[i] = moon_ecc
        if moon_semi < 0 or moon_ecc > 1 or moon_d < r_roche or moon_semi > 0.4*R_H:
          stable = False
      t += out_step
  except rebound.Collision as error:
    collision_happened = True
  end_wall_time = time.time()
  total_wall_time = end_wall_time - start_wall_time
  reason = "collision" if collision_happened else "ok"

  # Output data to file
  with open("output.txt", "a") as outfile:
    out_stg = "%1.2f, " % delta_H
    for i in range(0,n_moons):
        out_stg += "%1.4f, " % emax[i]
    out_stg += "%1.4e, %1.5f\n" % (sim.t/Pmoon, total_wall_time)
    outfile.write(out_stg)
  outfile.close()
  return

lock = mp.Lock()
phi = (1+np.sqrt(5))/2.
star_mass = 1 # Sun M
earth_mass = 3.003e-6 #Solar M
earth_semi = 1 #A_M
earth_ecc = 0.01671022 # Orbit Circulation
earth_radius = 6371e3/astronomical_unit
e_m = 1e-7
r_earth = 4.259e-5
k2_earth = 0.298
c_earth = 0.3308 # Dimensionless moment of inertia, normalized by MR^2
R_H = 1*(earth_mass/3)**(1/3.) #0.005

home = os.getcwd()
if moon_mass == "L":
    home += "/Luna/"
    M_moon = earth_mass/81. #Luna moon
    rho_moon = 3.3
    radius_moon = 1737e3/astronomical_unit
elif moon_mass == "C":
    home += "/Ceres/"
    M_moon = 4.72e-10 #Ceres mass
    rho_moon = 2.2
    radius_moon = 473e3/astronomical_unit
elif moon_mass == "P":
    home += "/Pluto/"
    M_moon = 0.002183*earth_mass #Pluto mass
    rho_moon = 1.86
    radius_moon = 1.19e6/astronomical_unit

r_roche = 2.44*r_earth*(5.515/rho_moon)**(1/3.) #Roche Radius 1.233e-04 #Luna 3.3 Ceres 2.2 Pluto 1.86
a_m = 2 * r_roche

#setup output file
header = "#delta_H, "
for i in range(n_moons):
    header += "max e_m" + str(i+1) + ", "
header += "sim finish time (n Pmoon), wall_time\n"
write_to_file("output.txt", header, "head")
for i in range(0, n_moons):
  fname = "evolution%i.txt" % i
  with open(fname, "w") as outfile:
    outfile.write("#Time, Semimajor, Ecc, Distance\n")
    outfile.close()
delta_rng = np.arange(8.59,8.60,0.1)
# delta_rng = np.concatenate((np.arange(3.5,4.5,0.01), np.arange(6.9,7.2,0.01)))
# for d in range(0,len(delta_rng)):
#   run_simulation(d)
pool = mp.Pool(processes=8)
pool.map(run_simulation, range(len(delta_rng)))
pool.close()
