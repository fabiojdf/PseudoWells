import matplotlib.pyplot as plt
from numpy import array, arange, linspace, full


def transition_matrix_show(matrix):
    fig, ax = plt.subplots(figsize=(5, 5))

    v = ax.matshow(array(matrix), cmap=plt.cm.Blues)
    for i in range(len(array(matrix))):
        for j in range(len(array(matrix))):
            c = array(matrix)[j, i]
            ax.text(i, j, str(c), va='center', ha='center')
    plt.colorbar(v, label='Probability')
    ax.set_title('Transition matrix', fontsize=12)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

    plt.savefig('results/transition_matrix.png', dpi=250, bbox_inches='tight')


def bivariate_plot(X, Y, Z, color):

    fig, ax = plt.subplots(figsize=(6, 5))

    for i in range(len(color)):
        cs = ax.contour(X, Y, Z[i], levels=6, colors=f'{color[i]}')

    ax.set_title('Bivariate Gaussian')
    ax.set_xlabel('PhiE [dec]')
    ax.set_ylabel('Vclay [dec]')
    ax.set_xticks(arange(0, 0.41, 0.1))
    ax.set_yticks(arange(0, 0.81, 0.2))
    ax.set_xlim(0., 0.4)
    ax.set_ylim(0, 0.8)
    ax.grid(which='major', color='#CCCCCC', linestyle='--', alpha=0.3)

    plt.savefig('results/bivariate_gaussian.png', dpi=250, bbox_inches='tight')


def variogram_plot(variance, range, sill, model="Exponential model"):
    
    fig, ax = plt.subplots(figsize=(6, 5))

    ax.plot(linspace(0, len(variance), len(variance)), variance, color='k', label=model)
    ax.axvline(x=range, color='blue', linestyle='--')
    ax.axhline(y=sill, color='blue', linestyle='--')
    ax.text(25, sill + sill/100, 'Sill')
    ax.text(range + range/100, 0.001, 'Range', rotation=-90)
    ax.set_title('Variogram')
    ax.set_xlabel('h [m]')
    ax.set_ylabel('Spatial variance (h)')
    ax.set_ylim(0, 0.0025)
    ax.set_xlim(0, 50)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.grid(which='major', color='#CCCCCC', linestyle='--', alpha=0.3)
    ax.legend(loc='lower right')

    plt.savefig('results/variogram.png', dpi=250, bbox_inches='tight')


def cov_matrix_plot(matrix):

    fig, ax = plt.subplots(figsize=(5, 5))

    v = ax.matshow(array(matrix), cmap=plt.cm.jet)
    plt.colorbar(v, label='Probability')
    ax.set_title('Covariance matrix', fontsize=12)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

    plt.savefig('results/covariance_matrix.png', dpi=250, bbox_inches='tight')


def rpm_plot(x0, y0, x1, y1, c0, c1, lith0, lith1):
    
    fig, ax = plt.subplots(figsize=(6, 5))

    ax.plot(x0, y0, color=c0, label=lith0)
    ax.plot(x1, y1, color=c1, label=lith1)
    ax.set_xlim(0, 0.4)
    ax.set_xlabel('PhiE [dec]')
    ax.set_label('K [GPa]')
    ax.set_title('Rock physics models')

    ax.legend(loc='upper right')

    plt.savefig('results/rpm.png', dpi=250, bbox_inches='tight')


def logplots(depth, vclay, phie, sw, Vpsat, Vssat, rhor, facies, facies_code, colors):

    fig, (ax, ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 7, sharey=True, figsize=(14, 10))

    for i in range(len(colors)):
        ax.fill_betweenx(depth, 0, 5.0, where=facies_code[i], color=colors[i], label=facies[i],interpolate=True)

    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax.set_title('Lithologies', fontsize=12)
    ax.set_xlim(0, 1.0)
    ax.legend(loc='lower left')

    ax1.plot(vclay, depth, color='green', lw=1.6)
    ax1.set_title('Vclay [dec]', fontsize=12)
    ax1.set_xticks([0.0, 0.5, 1.0])
    ax1.set_xlim(0, 1.0)

    ax2.plot(phie, depth, color='black', lw=1.6)
    ax2.set_title('PhiE [dec]', fontsize=12)
    ax2.set_xticks([0.0, 0.2, 0.4])
    ax2.set_xlim(0, 0.4)

    ax3.plot(sw, depth, '--k', lw=1.6)
    ax3.fill_betweenx(depth, 0, sw, color='deepskyblue', interpolate=True)
    ax3.set_title('SW [dec]', fontsize=12)
    ax3.set_xticks([0.0, 0.5, 1.0])
    ax3.set_xlim(0, 1.0)

    ax4.plot(Vpsat, depth, color='magenta', lw=1.6)
    ax4.set_title('Vp [km/s]', fontsize=12)
    ax4.set_xticks([1.5, 2.0, 3.0, 4.0, 5.0])
    ax4.set_xlim(1.5, 5.0)

    ax5.plot(Vssat, depth, color='purple', lw=1.6)
    ax5.set_title('Vs [km/s]', fontsize=12)
    ax5.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5])
    ax5.set_xlim(0.5, 2.5)

    ax6.plot(rhor, depth, color='blue', lw=1.6)
    ax6.set_title('RHOB [g/cm3]', fontsize=12)
    ax6.set_xticks([1.95, 2.5, 2.95])
    ax6.set_xlim(1.95, 2.95)

    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.grid(which='minor', alpha=0.2)

    plt.gca().invert_yaxis()
    plt.subplots_adjust(wspace=0.25)
    plt.savefig('results/logplot.png', dpi=250, bbox_inches='tight')


def crossplot_vpvs(vp, vs, z="blue", cmp=plt.cm.jet, zlabel=None):

    fig, ax = plt.subplots(figsize=(6,5))

    map = ax.scatter(vp, vs, c=z, s=6.5, cmap=cmp)
    ax.set_xlabel('Vp [km/s]', fontsize=12)
    ax.set_ylabel('Vs [km/s]', fontsize=12)
    ax.set_title('Crossplot Vp/Vs', fontsize=14)
    ax.set_xlim(2.0, 5.0)
    ax.set_ylim(1.0, 2.5)
    ax.set_xticks([2.0, 3.0, 4.0, 5.0])
    ax.set_yticks([1.0, 1.5, 2.0, 2.5])
    ax.grid(which='minor', alpha=0.4)
    
    if type(z) is not str:
        plt.colorbar(map, label=zlabel)

    plt.savefig('results/crossplot_vpvs.png', dpi=250, bbox_inches='tight')
