import TransitProject.webScraping as ws
import argparse

parser = argparse.ArgumentParser(description="Plot simulations and models")
parser.add_argument("--star", help="The name of the star to search the lightcurve of.", required=True)
parser.add_argument("--saveTTV", help="Whether the results of the TTV should be saved after graphing.", action="store_true")
args = parser.parse_args()

# if we are saving the TTV
if args.saveTTV:
    # Run the fetch & save code
    ws.fetchTESSTTV(target=args.star)
# otherwise, just plot
ws.plotTESS(target=args.star)
