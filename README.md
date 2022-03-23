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
│ ├ Poster.pdf
│ └ ProjectProposal.pdf
├ Source/
│ ├ old_code /
│ │ └ *Contains various bits legacy code that is kept for reference*
│ └ TransitProject/
└ README.md
```

In the `Source/` folder is all of the python scripts to be ran.
Additionally, contained within the `Source/TransitProject` folder is the individual code modules common across the project.

# Overview

The (Current) High-level plan for this project is (as everything else is) subject to change. It is as follows:

- Find suitable observation targets (See [Project Proposal](#project-proposal)).
- Observe targets to collect data.
- Fit data to transit model and upload to Exoclock (See [Literature Review](#literature-review)).
- Summarise methodology (See [Poster](#poster)).
- Write code.
  - Select target planetary system
  - Find all ExoClock/ETD data for the planetary system.
  - Find the known system information (ie: NASA Exoplanet Archive search).
  - Simulate the known parts of the system.
  - Fit the simulation to the ExoClock/ETD data.
  - Determine whether the system does exhibit TTV.
  - Run hundreds of simulations to suggest possible system layout.
  - Summarise every stage with graphs and charts.
- Compile results (see [Presentation](#presentation)).
- Summarise findings and results, and suggest options for future work (see [Dissertation](#dissertation)).

# Programming steps

This is a general overview of the "Write Code" step above, and is of course subject to change.

- Collect all required data
  - ETD (Lightcurves)
  - Exoclock (Lightcurves & system info)
  - Exoplanet Archive (System info)
  - Exoplanet Catalogue (System info & digital renders for pretty graphs)
  - others too
- Combine data
  - Merge ETD & Exoclock data so as to plot hisotrical lightcurves
  - Merger Exoplanet Archive & Exoplanet Catalogue & ExoClock to find system information
    - Figure out how to deal with conflicting information
- Simulation
  - Write a shared library for simulation
  - Setup a pipeline with multithreading instances
  - Allow the inclusion of other factors (Ie: maybe we want to look at the Effect of GR?)
  - Determine a standardised output
- Processing & Fitting
  - Determine a standardised input for data
  - Convert all data to standardised form & storage (ie: JSON)
  - Use multiple fitting types, ie: Idealised Sine, Sinusoid from known parameters, simple simulation
  - Fit using some extensive methodology, ie: MCMC, Simple Polyfit
- Analysis
  - Determine accuracy of fit to known data.
  - Probably run Fourier analysis.
  - Determine a way of pulling out key frequency terms.
  - Map key information to system parameters. ie: frequency to orbital period mapping.
- Re-simulate
  - Use analysis to construct a range of possible system layouts
  - Simulate each system layout (See Simulation step)
  - Re-perform analysis of simulation, determine which, if any, simulated system layouts fit the data.
- Pretty outputs
  - Graphical representations of all system layouts.
  - Annotated charts of any fit data.
  - Nice system parameter layout à la [wikipedia](https://en.wikipedia.org/wiki/Earth).
