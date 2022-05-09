# python modules
from matplotlib.lines import Line2D
import matplotlib.pylab as plt
import argparse, rebound

# my modules
import TransitProject.Simulation as ts
import TransitProject as tp

def fetchArgs():
    ''' fetch the program arguments from command line input. '''
    parser = argparse.ArgumentParser(description="Plot system configurations.")
    # select planet
    parser.add_argument("--planet", help="The name of the system to plot.", required=True)
    args = parser.parse_args()
    # return the argument input
    return args

# start the script
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()

    # fetch the parameters for the system
    df, params = ts.fetchParams(args.planet)

    # construct the simulation parameters
    simArray = ts.constructSimArray(df, params, useerror=False)

    # fetch the simulation for this
    sim = ts._arrayToSim_(simArray[0], params)

    # plot the simulation
    # plot the position of bodies
    fig, ax_main, ax_top, ax_right = rebound.OrbitPlot(sim, slices=0.5, xlim=[-2.,2], ylim=[-2.,2], color=True, unitlabel="[AU]")
    fig.suptitle(df.iloc[0]["name"])
    # create the legend
    cols = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]
    leg_elem = [Line2D([0], [0], color=cols[i], lw=1.5, label="{} Orbit".format(df.iloc[i+1]["name"])) \
                for i, item in enumerate(df.iloc[1:].index)]
    # apply the legend
    plt.legend(handles = leg_elem, framealpha=0, bbox_to_anchor=(0, 1.5), loc="upper left")
    # and save
    fig.tight_layout()
    fig.savefig("{}_layout.png".format(df.iloc[0]["name"]), transparent=False, bbox_inches='tight')
    fig.savefig("{}_transparent_layout.png".format(df.iloc[0]["name"]), transparent=True, bbox_inches='tight')
    plt.show()
