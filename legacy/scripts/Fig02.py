import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable

file_T0 = "/content/output_L_2moons_T0_9.0_0.txt"
file_T1 = "/content/output_L_2moons_T1_9.0_100.txt"
file_T2 = "/content/output_L_2moons_T1_9.0_698.txt"


earth_mass = 3.003e-6 #Solar M
M_moon = earth_mass/81 #L earth_mass/81, C 4.72e-10, P 0.002183*earth_mass

def MMR(y):
    X = 0.5 * ((2 * M_moon) / (3 * (earth_mass + M_moon)))**(1./3.)
    y = y**(2.0/3.0)
    beta = (1.0 / X) * (y - 1.0) / (y + 1.0)
    return beta

def plot_emax(ax, file, moon, idx):

    # Data from file
    data = np.genfromtxt(file, delimiter=",", comments="#",autostrip=True,
                         dtype=None,names=True,encoding='utf-8')
    DeltaH = data['dH']
    e0     = data['e0']
    fate   = data['fate']
    if moon == 1:
        Z_raw = data['eM1_max'] - data['eM1_min'] #M1 col 2
    elif moon == 2:
        Z_raw = data['eM2_max'] - data['eM2_min']#M2 col 5
    else:
        print("0")
    conditions = np.isin(fate, ['a<0','a>0.4R_H','collision','d<R_roche'])
    Z_raw[conditions] = 1.5

    fs = 20
    lbl_fs = 40
    lw = 8
    sublbl = ['a','b','c','d','e','f']
    #data grid points
    xi = np.arange(3.5, 9.01, 0.01) #x axis
    yi = np.arange(0, 0.31, 0.01) # y axis
    Xi, Yi = np.meshgrid(xi, yi)

    Zi = griddata((DeltaH, e0), Z_raw, (Xi, Yi), method='nearest')
    Zi = np.ma.masked_invalid(Zi)

    #color map
    cmap = cm.gnuplot
    vmin, vmax = -2, 0

    norm = colors.Normalize(vmin, vmax)
    cmmapable = cm.ScalarMappable(norm, cmap)
    cmmapable.set_array([])

    ax.pcolormesh(xi, yi, np.log10(Zi), cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    #moon idx [0,1,2,3]
    if idx in [4, 5]:
        ax.set_xlabel("$\\beta \ (R_{H,m})$", fontsize=lbl_fs)
    else:
        ax.set_xticklabels([])
    if idx == 0 or idx == 2 or idx == 4:
        ax.set_ylabel("$e_0$", fontsize=lbl_fs)
    else:
        ax.set_yticklabels([])
    #tick params, def
    ax.minorticks_on()
    ax.tick_params(which='major', axis='both', direction='out', length=12.0, width=4.0, labelsize=fs)
    ax.tick_params(which='minor', axis='both', direction='out', length=6.0, width=2.0)
    ax.text(0.025, 0.92, sublbl[idx], color='k', fontsize='xx-large', weight='bold',
            horizontalalignment='left', transform=ax.transAxes)
    ax.set_xlim(3.5, 9.0)
    ax.set_ylim(0, 0.3)
    if idx == 5:
        color_label = "$\log_{10}\Delta e$"
        cax = fig.add_axes([0.92, 0.11, 0.015, 0.77])
        cbar = plt.colorbar(cmmapable, cax=cax, orientation='vertical')
        cbar.set_label(color_label, fontsize=lbl_fs)
        cbar.ax.tick_params(axis='both', direction='out', length=8.0, width=4.0, labelsize=fs)
    top_ticks =[5.9678,5.3222,6.2119,5.6769,4.8769,4.2948,4.6070,3.9278,8.75]
    top_ticks_lbl = ['8:1','6:1','9:1','7:1','5:1','4:1','9:2','7:2','$\\beta_\max$']
    if idx in [0,1]:
        ax_top = ax.twiny()
        ax_top.set_xlim(3.5, 9.0)
        ax_top.set_xticks(top_ticks)
        ax_top.set_xticklabels([])
        ax_top.tick_params(which='major', axis='both', direction='out', length=12.0, width=4.0, labelsize=fs)
        ax_top.tick_params(which='minor', axis='both', direction='out', length=6.0, width=4.0)
    for label, x in zip(top_ticks_lbl, top_ticks):
        ax.axvline(x=x, color='k', linestyle='--', linewidth=2)
        if idx in [0,1]:
            ax_top.set_xticklabels(top_ticks_lbl, fontsize=10,ha='center',
                    va='bottom', rotation = 45)


width = 10.
aspect = 2.1
fig = plt.figure(figsize=(14, 12), dpi=300)
#ax def
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)
# def moon, indx , files etc
#0 row 1
plot_emax(ax1, file_T0, moon=1, idx=0)
plot_emax(ax2, file_T0, moon=2, idx=1)

#100 row 2
plot_emax(ax3, file_T1, moon=1, idx=2)
plot_emax(ax4, file_T1, moon=2, idx=3)

#698 row 3
plot_emax(ax5, file_T2, moon=1, idx=4)
plot_emax(ax6, file_T2, moon=2, idx=5)

fig.subplots_adjust(wspace=0.1, hspace=0.15, right=0.90)
fig.savefig("L-emax.png", bbox_inches='tight', dpi=300)