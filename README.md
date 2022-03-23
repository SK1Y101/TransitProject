# Project
A repository for my Masters Project on exoplanet transit timing variations.

[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)

[![Website Link](https://img.shields.io/badge/Website-Link-aqua?labelColor=lightblue&style=for-the-badge)](https://sk1y101.github.io/projects/TransitProject/)


![GitHub](https://img.shields.io/github/license/SK1Y101/TransitProject)
[![CodeFactor](https://www.codefactor.io/repository/github/SK1Y101/TransitProject/badge)](https://www.codefactor.io/repository/github/SK1Y101/TransitProject)
[![wakatime](https://wakatime.com/badge/github/SK1Y101/TransitProject.svg)](https://wakatime.com/badge/github/SK1Y101/TransitProject)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/SK1Y101/TransitProject)
![GitHub last commit](https://img.shields.io/github/last-commit/SK1Y101/TransitProject)

![GitHub language count](https://img.shields.io/github/languages/count/SK1Y101/TransitProject)
![GitHub top language](https://img.shields.io/github/languages/top/SK1Y101/TransitProject)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/SK1Y101/TransitProject)
![Lines of code](https://img.shields.io/tokei/lines/github.com/SK1Y101/TransitProject)

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
│ ├ TransitProject/
│ │ ├ __init__.py
│ │ └ webScraping.py
│ └ fetchData.py
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

This is a general overview of the "Write Code" step above, and is of course subject to change. Next to each item is a tickbox that marks its completion status.

- [ ] Collect all required data
  - [ ] ETD Lightcurves
  - [ ] Exoclock
    - [ ] Lightcurves or mid-transit times
    - [ ] system info
  - [ ] Exoplanet Archive System info
  - [ ] Exoplanet Catalogue
    - [ ] System info
    - [x] digital renders for pretty graphs (Completed from another of my projects, available [here](https://github.com/SK1Y101/Data_Collection_Pipeline))
  - [ ] other data sources
- [ ] Combine data
  - [ ] Merge ETD & Exoclock data so as to plot hisotrical lightcurves
  - [ ] Merger Exoplanet Archive & Exoplanet Catalogue & ExoClock to find system information
    - [ ] Figure out how to deal with conflicting information
- [ ] Simulation
  - [ ] Write a shared library for simulation
  - [ ] Setup a pipeline with multithreading instances
  - [ ] Allow the inclusion of other factors (Ie: maybe we want to look at the Effect of GR?)
  - [ ] Determine a standardised output
- [ ] Processing & Fitting
  - [ ] Determine a standardised input for data
  - [ ] Convert all data to standardised form & storage (ie: JSON)
  - [ ] Use multiple fitting types, ie: Idealised Sine, Sinusoid from known parameters, simple simulation
  - [ ] Fit using some extensive methodology, ie: MCMC, Simple Polyfit
- [ ] Analysis
  - [ ] Determine accuracy of fit to known data.
  - [ ] Probably run Fourier analysis.
  - [ ] Determine a way of pulling out key frequency terms.
  - [ ] Map key information to system parameters. ie: frequency to orbital period mapping.
- [ ] Re-simulate
  - [ ] Use analysis to construct a range of possible system layouts
  - [ ] Simulate each system layout (See Simulation step)
  - [ ] Re-perform analysis of simulation, determine which, if any, simulated system layouts fit the data.
- [ ] Pretty outputs
  - [ ] Graphical representations of all system layouts.
  - [ ] Annotated charts of any fit data.
  - [ ] Nice system parameter layout à la [wikipedia](https://en.wikipedia.org/wiki/Earth).
