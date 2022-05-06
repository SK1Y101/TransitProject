import TransitProject.webScraping as ws
import argparse

parser = argparse.ArgumentParser(description="Plot simulations and models")
parser.add_argument("--star", help="The name of the star to search the lightcurve of.", required=True)
args = parser.parse_args()

ws.plotTESS(target=args.star)
