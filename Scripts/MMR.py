earth_mass = 3.003e-6 #Solar M
M_moon = earth_mass/81 #L earth_mass/81, C 4.72e-10, P 0.002183*earth_mass

def MMR(y):
    X = 0.5 * ((2 * M_moon) / (3 * (earth_mass + M_moon)))**(1./3.)
    y = y**(2.0/3.0)
    beta = (1.0 / X) * (y - 1.0) / (y + 1.0)
    return beta
# for i in range()
print(MMR(24/1))