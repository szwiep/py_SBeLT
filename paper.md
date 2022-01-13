---
title: 'pyBeLT: A Python software for fluvial sediment transport under rarefied conditions'
tags:
  - Python
  - geomorphology
  - sediment transport
  - stochastic
  - Poisson
authors:
  - name: Sarah Zwiep^[co-first author] 
    orcid: 
    affiliation: "1"
  - name: Shawn M. Chartrand^[co-first author] 
  - orcid: 0000-0002-9309-1137
    affiliation: 1
affiliations:
 - name: School of Environmental Science, Simon Fraser University
   index: 1
date: 12 January 2022
bibliography: paper.bib

---

# Summary

Granular sediment of various sizes moves downstream along river beds when water flow is capable of entraining particles from 
the bed surface. This process is known as bed load sediment transport because the particles travel close to the boundary. 
It is common to treat the transport process as a predictive problem in which the mean transport rate past a stationary 
observation point is a function of local water flow conditions (Parker, 2008; Wainwright, 2014; Ancey, 2020). However, 
deterministic approaches to the problem neglect the stochastic nature of transport, which originates from the movements 
of individual particles (Einstein, 1937; Furbish and Doane, 2021). Here, we present an open-source Python software, pySBeLT, 
which simulates the kinematics of rarefied particle transport (low rates) as a stochastic process along a riverbed profile. 
This software is motivated by a need to better understand connections between individual particle motions and local transport 
rates, or the flux.

# Statement of need

Research at the intersection of geomorphology, geophysics and hydraulics is increasingly focused on building 
a theoretical foundation for the treatment of bedload transport as a stochastic phenomenon (e.g. Ancey, 2020a; 
Furbish and Doane, 2021). Associated theories are commonly tested against laboratory data from "rarefied" transport 
conditions (Furbish et al., 2016), where transport rates are low to moderate, interactions between two or more 
moving particles are rare, and a relatively small fraction of particles on the bed surface participate in transport 
(e.g. Ancey, 2010; Roseberry et al., 2012; Fathel et al., 2015; Wu et al., 2019). For example, laboratory experiments 
using a downstream light table counting device and conducted at roughly twice the threshold for particle motion involve 
the transport of less than approximately 12% of particles on the upstream bed surface (Chartrand, 2017). This result 
highlights that the flux measured across a boundary or within an area of bed surface is directly linked to the motions 
of individual particles arriving from upstream locations (Furbish et al., 2012). 

Because particle movements are controlled by fluid turbulence, the irregular bed surface, and collective movement effects 
(Ancey et al., 2006; Ancey, 2008; Lee and Jerolmack, 2019), the connection between particle movements and the bedload 
transport rate has been difficult to formulate mathematically. The pySBeLT software provides an extensible framework within 
Python to numerically examine correlations between upstream particle entrainment rates and travel distances, with downstream 
flux. The pySBeLT software was motivated by a birth-death, immigration-emigration Markov model for bedload transport. Here, 
the movements of individual particles are represented by stochastic entrainment, motion, and deposition processes, and sediment 
flux is represented as a counting phenomenon where the number of particles in motion above the bed surface is a random 
variable (Ancey 2008). The pySBeLT software supports ensemble simulations so that repeat numerical experiments can be conducted,
or the problem can be efficiently probed across a range of input parameter values (discussed below).

pyBeLT is run forward in time according to default or user specified parameter values in param.yaml (see the README.md for 
more details). After initialization, pyBeLT first constructs a bed of fixed particles of set_diam in both the downstream and 
cross-stream dimensions (one particle wide in the present build), and over a downstream domain length x_max. Bed surface particles
of set_diam are then randomly placed at vertices between fixed bed particles until the pack_density is met. Vertices are defined 
by a contact point between two adjacent particles. The bed of surface particles is then separated into num_subregions, and at this 
point the forward simulations are ready to commence. Simulation iterations involve three steps: (1) the number of particle entrainment 
events per num_subregions are drawn from a Poisson pmf, and this is done randomly for each numerical step up to n_iterations; 
(2) surface particles from each subregion are randomly selected for entrainment, and if there are insufficient surface particles available 
for entrainment, then all available particles are moved; (3) each entrained particle moves a distance according to a randomly sampled value 
from either the normal or lognormal distribution, and is placed at the nearest vertex between two particles that is available for placement. 
Placed particles are permitted to stack up to the level_limit in height. Travel distances of particles that exceed x_max are returned and 
queued at the upstream boundary, and are introduced back into the domain at the next numerical step according to travel distance sampling 
described above. This overall process repeats for the specified n_iterations. 

pyBeLT tracks a number of different parameters through a simulation: the vertical and horizontal positions of every particle center, 
the randomly sampled number of entrainment events, the number of particles actually entrained, the randomly sampled particle travel 
distance, the actual particle travel distance, the particle ‘age’, or the number of numerical steps since last entrainment for every 
particle, and the number of particles which cross all boundaries, i.e. sub-region and downstream at x_max. All values are written to ____ files. 
pyBeLT produces a time varying signal of particle flux counted at the downstream domain (as well as internal subregion domains), with a 
particle bed that changes through particle stacking and pile removal, and downstream motions of travel distance (Fig. 1). An implication 
of particle stacking within the context of the pySBeLT stochastic framework is a time varying signal of the average “particle age”, as well 
as the average “particle age range”, defined as the difference of the maximum and minimum particle ages. The software can be readily modified 
to simulate kinematics using different probability distributions, or examining particle age dynamics for deeper beds of particles available 
for transport. py_BeLT can also be extended to 2-dimensions in the cross-stream. The relatively simple parameterization of pySBeLT execution 
also makes it suitable for use as a teaching tool within advanced undergraduate and graduate courses emphasizing bedload transport.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures


|![Image](/figures/Figure%201.png)
|:--:| 
| *Example py_BeLT output of particle flux at downstream boundary and particle bed configuration at numerical step 100* |

# Acknowledgements

S. Zwiep was funded in part through an Undergraduate Student Research Award from the 
National Science and Engineering Research Council of Canada (NSERC). S.M. Chartrand was 
funded through a Postdoctoral Fellowship awarded by NSERC, and through internal research 
funding provided by Simon Fraser University. The model was inspired by discussions with 
David Jon Furbish, who also provided useful input and critical feedback at various stages 
of model development and testing. Kevin Pierce also provided helpful feedback during model 
development. 


# References
