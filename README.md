# Project
A repository for my Masters Project on exoplanet transit timing variations.

[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)

[![Website Link](https://img.shields.io/badge/Website-Link-aqua?labelColor=lightblue&style=for-the-badge)](https://sk1y101.github.io/projects/TransitProject/)

# Motivation
Transiting exoplanets are growing ever more numerous, since *Kepler* launched in 2009, we have seen an explosive growth in the number of them detected. There are thousands visible in any single patch of the sky, and that's awesome.

What interests me, however, is the potential hidden planets in those systems. The planets whose transit is off axis to our line of sight, and slip past silent in the night. Through their gravitational influence on the planets we can see, hints of their presence are evident, hints that we will attempt to locate programmatically.

# File Structure

The file structure of this repository is as follows:

```
Repo/
├ Files/
│ ├ Literature_Review.pdf
│ └ ProjectProposal.pdf
├ Source/
│ ├ ProjectModules/
│ │ ├ trig_funcs.py
│ │ └ nbody.py
│ └ multi-planetary_transit.py
└ README.md
```

In the `Source/` folder is all of the python scripts to be ran.
Additionally, contained within the `Source/ProjectModules` folder is the individual code modules common across the project.

ie: `nbody.py` takes the equations of motion and allows additional code to solve them, while `trig_funcs.py` is a set of functions for handling trigonometry to arbitrary levels of precision (Using Decimal)
