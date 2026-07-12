import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import astronomical_unit
from matplotlib.ticker import FixedLocator, NullFormatter
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

lbl_fs = 'x-large'
colors = ["r","b"]
labels = ["m1", "m2"]
tau = [0, 100, 698]
sublbl = ['a','b', 'c', 'd', 'e','f']
matplotlib.rcParams["pdf.compression"] = 9
matplotlib.rcParams["savefig.dpi"] = 300
fig, axes = plt.subplots(2, 3, figsize=(14, 6), dpi=300)

for col, i in enumerate(tau):

    data = np.loadtxt("Luna_%i_final.txt" % i, delimiter=",", skiprows=1, usecols=(0, 2, 5, 8))
    delta_H = data[:, 0]
    time = data[:, 3]

    axes[0, col].scatter(delta_H, time, s=1, color="black")
    axes[0, col].set_yscale("log")
    axes[0, col].axvline(x=8.75, color='k', linestyle='dashed', linewidth=1)
    axes[0, col].set_xlim(2 * np.sqrt(3), 9)
    axes[0, col].set_ylim(1e3, 1.5e7)
    axes[0, col].set_ylabel("$\\log_{10}{t}$" if col == 0 else "", fontsize=lbl_fs)
    axes[0, col].minorticks_on()
    axes[0, col].tick_params(which='both', direction='out', length=6, width=3, labelsize=12)
    axes[0, col].tick_params(which='minor', direction='out', length=3, width=2)
    axes[0, col].text(0.05, 0.9, sublbl[col], color='k', fontsize=15,
                      weight='bold',horizontalalignment='left', transform=axes[0, col].transAxes)
    if col != 0:
        axes[0, col].tick_params(labelleft=False)
        axes[1, col].tick_params(labelleft=False)
    axes[0, col].tick_params(labelbottom=False)

    for j in range(1, 3):
        e = data[:, j]
        axes[1, col].scatter(delta_H, e,
            s=2,color=colors[j-1],label=labels[j-1])
    axes[1, col].set_yscale("log")
    axes[1, col].set_xlim(2 * np.sqrt(3), 9)
    axes[1, col].set_ylim(1e-3, 1.3)
    axes[1, col].set_ylabel("$e_{max}$" if col == 0 else "",fontsize=lbl_fs)
    axes[1, col].minorticks_on()
    axes[1, col].tick_params(which='both', direction='out', length=6, width=3, labelsize=12)
    axes[1, col].tick_params(which='minor', direction='out', length=3, width=2)
    axes[1, col].text(0.05, 0.9, sublbl[3 + col], color='k', fontsize=15,
                    weight='bold',horizontalalignment='left', transform=axes[1, col].transAxes)
    axes_top = axes[0, col].twiny()
    axes_top.set_xlim(axes[0, col].get_xlim())
    axes_top.set_xticks([8.75])
    axes_top.set_xticklabels(['$\\beta_{max}$'], fontsize=10,ha='center',
                    va='bottom', rotation = 45)
    axes[1, col].set_xlabel("$\\beta \ (R_{H,m})$", fontsize=lbl_fs)

    if col == 1: #and col == 2 and r == 2:
        axins = axes[1, col].inset_axes([0.11, 0.2, 0.35, 0.35]) #left bottom width heght
        for j in range(1, 3):
            e = data[:, j]
            axins.scatter(delta_H, e,
                s=2,color=colors[j-1],label=labels[j-1])
        # axins[1, col].scatter(delta_H, time, s=1, color='r', alpha=0.7)
        axins.set_yscale('log')   # zoomed semimajor-axis range
        axins.set_xlim(7, 7.6)          # zoomed time range
        axins.set_ylim(1e-1, 1.8e-1 ) 
        mark_inset(axes[1, col], axins,loc1=2, loc2=4,
                    fc="none", ec="k", lw=1)
        # axins.tick_params(labelsize=8)
        axins.yaxis.set_major_locator(FixedLocator([1e-1]))
        axins.set_yticklabels(["$10^{-1}$"])
        axins.yaxis.set_minor_locator(FixedLocator([1.1e-1, 1.2e-1, 1.3e-1, 1.4e-1,
                                                    1.5e-1, 1.6e-1, 1.7e-1, 1.8e-1]))
        axins.yaxis.set_minor_formatter(NullFormatter())
        axins.tick_params(axis='y', which='major', length=5, width=1.2, labelsize=7)
        axins.tick_params(axis='y', which='minor', length=2.5, width=0.8)
        axins.tick_params(axis='x', which='both', labelsize=7)

#Color Bar
cmap = mpl.colors.ListedColormap(colors)

bounds = np.arange(len(colors) + 1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cax = fig.add_axes([0.125, 0.92, 0.08, 0.03])

cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                  cax=cax, orientation='horizontal', drawedges=True)

cax.xaxis.set_ticks_position('top')
cb.set_ticks(bounds[:-1] + 0.5)
cb.set_ticklabels([str(i+1) for i in range(len(colors))], fontsize=12)
cax.tick_params(axis='x', length=0)

cax.text(1.02, 0.5, "Moon Index $i$", transform=cax.transAxes,
         va='center', ha='left', fontsize=14)

fig.subplots_adjust(wspace=0.1, hspace=0.15)#, right=0.90, top=0.86)
plt.savefig("logt_ecc_Luna", dpi=300, bbox_inches="tight", pad_inches=0.05)
plt.show()
