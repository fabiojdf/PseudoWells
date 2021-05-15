import matplotlib.pyplot as plt
from numpy import array

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


def logplots(depth, vclay, phie, sw, Vpsat, Vssat, rhor, shale, silty, sand):

    fig, (ax, ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 7, sharey=True, figsize=(14, 10))

    ax.fill_betweenx(depth, 0, 5.0, where=shale, color='green', label='Shale',interpolate=True)
    ax.fill_betweenx(depth, 0, 5.0, where=silty, color='brown', label='Silty-sand',interpolate=True)
    ax.fill_betweenx(depth, 0, 5.0, where=sand, color='yellow', label='Sand',interpolate=True)
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
    ax4.set_xticks([2.0, 4.0, 6.0])
    ax4.set_xlim(2, 6.0)

    ax5.plot(Vssat, depth, color='purple', lw=1.6)
    ax5.set_title('Vs [km/s]', fontsize=12)
    ax5.set_xticks([1.0, 2.5, 4.0])
    ax5.set_xlim(1.0, 4.0)

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
    ax.set_xlim(2.0, 6.0)
    ax.set_ylim(1.0, 4.0)
    ax.set_xticks([2.0, 4.0, 6.0])
    ax.set_yticks([1.0, 2.5, 4.0])
    ax.grid(which='minor', alpha=0.4)
    
    if type(z) is not str:
        plt.colorbar(map, label=zlabel)

    plt.savefig('results/crossplot_vpvs.png', dpi=250, bbox_inches='tight')
