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
<img src="https://www.openhub.net/p/TransitProject/widgets/project_thin_badge?format=gif" alt="https://www.openhub.net/p/TransitProject/widgets/project_thin_badge?format=gif" style="border-radius: 0.25rem;">

# Motivation
Transiting exoplanets are growing ever more numerous, since *Kepler* launched in 2009, we have seen an explosive growth in the number of them detected. There are thousands visible in any single patch of the sky, and that's awesome.

What interests me, however, is the potential hidden planets in those systems. The planets whose transit is off axis to our line of sight, and slip past silent in the night. Through their gravitational influence on the planets we can see, hints of their presence are evident, hints that we will attempt to locate programmatically.

# File Structure

The file structure of this repository is as follows:

```
Repo/
├ Files/
│ ├ Documents /
│ │ ├ Literature_Review.pdf
│ │ ├ Poster.pdf
│ │ ├ ProjectPresentation.pptx
│ │ └ ProjectProposal.pdf
│ └ Images /                    *Contains various images generated for use in the research paper*
├ Source/
│ ├ old_code /                  *Contains various bits legacy code that is kept for reference*
│ ├ raw_data/
│ │ ├ midTransitTimes/          *Contains CSV Data for all exoplanets mid-transit times*
│ │ ├ tess_data/                *Contains CSV Data and Light-curve fits for all exoplanet TESS Data*
│ │ ├ ETDEphemerides.csv        *ephemerides from Exoplanet Transit database*
│ │ ├ exoClockEphemerides.csv   *ephemerides from ExoClock database*
│ │ ├ exoplanetList.csv         *list of all exoplanets with mid-transit times*
│ │ ├ ps.csv                    *NASA Exoplanet Archive planetary information table*
│ │ └ pscomparrs.csv            *NASA Exoplanet Archive composite planetary table*
│ ├ sim_data/                   *Contains the outputs of any simulations that were run*
│ ├ TransitProject/
│ │ ├ __init__.py
│ │ ├ modelling.py
│ │ ├ Simulation.py
│ │ ├ tesslc.py
│ │ └ webScraping.py
│ ├ ControlCandidates.csv
│ ├ determineTTVCandidates.py
│ ├ fetchData.py
│ ├ processing.py
│ └ simulationPipeline.py
├ .gitignore
├ README.md
└ requirements.txt
```

In the `Source/` folder is all of the python scripts to be ran.
These scripts are:
 - `determineTTVCandidates.py` - Determine which, of all exoplanets, are ideal candidates for testing TTV Models.
 - `fetchData.py` - Fetch Ephemerides and mid-transit times for exoplanets.
 - `fetchFittingParams.py` - Generates Latex tables from the model fitting parameters.
 - `plotOrbitalConfig.py` - Plot the layout for any given exoplanetary system
 - `processing.py` - Fit a transiting model to the observed data.
 - `simulationPipeline.py` - Simulate the known properties of an exoplanetary system to fetch transit times.

Additionally, there are the following files in `Source/`:
 - `ControlCandidates.csv` - Comma separated list of ideal testing candidates, in decreasing order of mid-transit data-points

Contained within `Source/TransitProject/` is the individual code modules common across the project.
Additionally, within `Source/old_code/` is legacy code that may still be useful.

At some point, I will probably have to write some unit-tests for all of my methods to ensure we have the expected outputs. Unfortunately, although I know how useful those are, I don't like doing that, so that may be a while (or possibly never, I'm not sure)

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

- [x] Collect all required data
  - [x] ETD Lightcurves or mid-transit times
  - [x] Exoclock
    - [x] Lightcurves or mid-transit times
    - [x] ephemerides
  - [x] Exoplanet Archive System info
  - [x] Exoplanet Catalogue
    - [x] System info (Completed from another of my projects, available [here](https://github.com/SK1Y101/Data_Collection_Pipeline))
    - [x] digital renders for pretty graphs (Completed from another of my projects, available [here](https://github.com/SK1Y101/Data_Collection_Pipeline))
  - [x] TESS/Kepler raw data > fit lightcurves to that for comparisons.
  - [x] Refactor code
- [x] Combine data
  - [x] Merge ETD & Exoclock data so as to plot hisotrical lightcurves
  - [x] Merge Exoplanet Archive & Exoplanet Catalogue & ExoClock to find system information
    - [x] Figure out how to deal with conflicting information
  - [x] Refactor code
- [x] Simulation
  - [x] Write a shared library for simulation
    - [x] Simulation pipeline
      - [x] Fetch the planetary system to simulate
      - [x] Fetch the system parameters
      - [x] Simulate the system with all known bodies to determine perturbed transit times
      - [x] Simulate the above with all known errors to determine error bars on perturbed transit times
  - [x] Setup a pipeline with multithreading instances for many simulations in parallel
  - [x] Allow the inclusion of other factors (Ie: maybe we want to look at the Effect of GR?)
  - [x] Determine a standardised output
  - [x] Refactor code
- [x] Processing & Fitting
  - [x] Determine a simpler model that can produce the sinusoid for a given system layout
    - [x] Ensure this is accurate to system simulation
    - [x] Should be decently fast to run
    - [x] Is possible to give to an MCMC fitter
  - [x] Linear fit to the visual data
    - [x] Fit a line to the data
    - [x] compare against exoclock ephemerides
    - [x] Search through Residuals for TTV
  - [x] Use multiple fitting types, ie: Idealised Sine, Sinusoid from known parameters, simple simulation
  - [x] Fit using some extensive methodology, ie: MCMC, Simple Polyfit
  - [x] Refactor code
- [x] Analysis
  - [x] Determine accuracy of fit to known data.
    - [x] BIC Fit from bayesian analysis
  - [x] Run MCMC fitting using TTVFast
  - [x] Probably run Fourier analysis.
  - [x] Determine a way of pulling out key frequency terms.
  - [x] Map key information to system parameters. ie: frequency to orbital period mapping.
  - [x] Determine if the known parameters of the system adequately map to the observed variation
    - [x] If not, determine the parameters of the perturbation.
  - [x] Refactor code
- [x] Re-simulate
  - [x] Use analysis to construct a range of possible system layouts (See Analysis)
  - [x] Simulate each system layout (See Simulation step)
  - [x] Re-perform analysis of simulation, determine which, if any, simulated system layouts fit the data.
  - [x] Refactor code
- [x] Pretty outputs
  - [x] Graphical representations of all system layouts.
  - [x] Annotated charts of any fit data.
  - [x] Nice system parameter layout à la [wikipedia](https://en.wikipedia.org/wiki/Earth).
  - [x] Refactor code
- [x] Code Objectives Complete
