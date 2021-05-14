import numpy as np
import matplotlib.pyplot as plt


def logplots(depth, vclay, phie, sw, Vpsat, Vssat, rhor, shale, silty, sand):

    fig, (ax, ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 7, sharey=True, figsize=(10, 16))

    ax.fill_betweenx(depth, 0, 5.0, where=shale, color='green', interpolate=True)
    ax.fill_betweenx(depth, 0, 5.0, where=silty, color='brown', interpolate=True)
    ax.fill_betweenx(depth, 0, 5.0, where=sand, color='yellow', interpolate=True)
    ax.set_xlim(0, 1.0)

    ax1.plot(vclay, depth)
    ax1.set_xlim(0, 1.0)

    ax2.plot(phie, depth)
    ax2.set_xlim(0, 1.0)

    ax3.plot(sw, depth)
    ax3.set_xlim(0, 1.0)

    ax4.plot(Vpsat, depth)
    ax4.set_xlim(0, 4.0)

    ax5.plot(Vssat, depth)
    ax5.set_xlim(0, 4.0)

    ax6.plot(rhor, depth)
    ax6.set_xlim(1.95, 2.95)

    plt.gca().invert_yaxis()


    plt.tight_layout()

    plt.show()

    return fig