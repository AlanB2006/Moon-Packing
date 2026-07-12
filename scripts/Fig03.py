fig, axes = plt.subplots(3, 1, figsize=(4, 8), dpi=300)
col = [1, 2, 5]
s=2
R_alpha = 1
S_alpha = 0.5
reb_data = np.genfromtxt("reb_tide.txt",delimiter=',')
sec_data = np.genfromtxt("sec_tide.txt",delimiter=',')

for col, i in enumerate(col):
    ax = axes[col]
    if i == 5:
        reb_spin = reb_data[1:,5] / (2*np.pi*365.25)
        sec_spin = sec_data[1:,5] / (2*np.pi*365.25)

        ax.semilogx(reb_data[1:,0], reb_spin, 'r.', ms=s, alpha=R_alpha)
        ax.semilogx(sec_data[1:,0], sec_spin, 'b.', ms=s, alpha=S_alpha )
        ax.set_ylabel("$\Omega_{\mathrm {spin}}\;(\mathrm{day^{-1}})$")
        ax.set_xlabel("Time (yr)")
        ax.tick_params(labelbottom=True)
    else:
        ax.semilogx(reb_data[1:,0], reb_data[1:,i], 'r.', ms=s, alpha=R_alpha)
        ax.semilogx(sec_data[1:,0], sec_data[1:,i], 'b.', ms=s, alpha=S_alpha)
        ax.tick_params(labelbottom=False)

    if i == 1:
        ax.legend(["Reboundx", "Secular"])
        ax.set_ylabel("$a(R_{H})$")
    elif i == 2:
        ax.set_ylabel("$e_{sat}$")
    ax.minorticks_on()
    ax.tick_params(which='both', direction='out', length=6, width=3, labelsize=12)
    ax.tick_params(which='minor', direction='out', length=3, width=2)

fig.subplots_adjust(wspace=0.15, hspace=0.15)
plt.savefig("Reb_vs_Sec", dpi=300, bbox_inches="tight")