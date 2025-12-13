import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable

file_T0 = "output_L_2moons_T0.txt"
file_T1 = "output_L_2moons_T1.txt"

def plot_emax(ax, file, moon, idx):

    # Data from file
    data = np.genfromtxt(file, delimiter=",", comments="#")
    DeltaH = data[:, 0]
    e0     = data[:, 1]

    if moon == 1:
        Z_raw = data[:, 2] #M1 col 2
    elif moon == 2:
        Z_raw = data[:, 5] #M2 col 5
    else:
        print("0")

    fs = 20         
    lbl_fs = 40     
    lw = 8
    sublbl = ['a','b','c','d']

    xi = np.arange(3.5, 6.01, 0.01) #x axis
    yi = np.arange(0, 0.31, 0.01) # y axis
    Xi, Yi = np.meshgrid(xi, yi)
    
    Zi = griddata((DeltaH, e0), Z_raw, (Xi, Yi), method='linear')
    Zi = np.ma.masked_invalid(Zi)

    #color map 
    cmap = cm.gnuplot
    vmin, vmax = 0.001, 1.0
    
    norm = colors.Normalize(vmin, vmax)
    cmmapable = cm.ScalarMappable(norm, cmap)
    cmmapable.set_array([])
    
    ax.pcolormesh(xi, yi, Zi, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    #moon idx [1,2,3,4]
    if idx > 2:
        ax.set_xlabel("initial $\Delta_H$", fontsize=lbl_fs)
        
    if idx == 1 or idx == 3:
        ax.set_ylabel("initial $e_0$", fontsize=lbl_fs)

    ax.minorticks_on()
    ax.tick_params(which='major', axis='both', direction='out', length=12.0, width=4.0, labelsize=fs)
    ax.tick_params(which='minor', axis='both', direction='out', length=6.0, width=2.0)
    
    ax.text(0.90, 0.92, sublbl[idx-1], color='k', fontsize='xx-large', weight='bold', 
            horizontalalignment='left', transform=ax.transAxes)

    ax.set_xlim(3.5, 6.0)
    ax.set_ylim(0, 0.3)
    
    if idx == 4:
        color_label = '$max_e$'
        cax = fig.add_axes([0.92, 0.11, 0.015, 0.77])
        cbar = plt.colorbar(cmmapable, cax=cax, orientation='vertical')
        cbar.set_label(color_label, fontsize=lbl_fs)
        cbar.set_ticks(np.arange(0., 1.25, 0.25))
        cbar.ax.tick_params(axis='both', direction='out', length=8.0, width=2.0, labelsize=fs)


width = 10.
aspect = 2.1
fig = plt.figure(figsize=(14, 12), dpi=300)

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

plot_emax(ax1, file_T0, moon=1, idx=1)
plot_emax(ax2, file_T1, moon=1, idx=2)
plot_emax(ax3, file_T0, moon=2, idx=3)
plot_emax(ax4, file_T1, moon=2, idx=4)

fig.subplots_adjust(wspace=0.25, hspace=0.15, right=0.90)

fig.savefig("emax.png", bbox_inches='tight', dpi=300)
