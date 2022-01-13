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
observation point is a function of local water flow conditions [@Parker2008; @Wainwright2014; @Ancey2020]. However, 
deterministic approaches to the problem neglect the stochastic nature of transport, which originates from the movements 
of individual particles [@Einstein1937; @FurbDoane2021]. Here, we present an open-source Python model, pySBeLT, 
which simulates the kinematics of rarefied particle transport (low rates) as a stochastic process along a riverbed profile. 
This model is motivated by a need to better understand connections between individual particle motions and local transport 
rates, or the flux.

# Statement of need

Research at the intersection of geomorphology, geophysics and hydraulics is increasingly focused on building 
a theoretical foundation for the treatment of bedload transport as a stochastic phenomenon [@Ancey2020; 
@FurbDoane2021]. Associated theories are commonly tested against laboratory data from "rarefied" transport 
conditions [@Furb2016], where transport rates are low to moderate, interactions between two or more 
moving particles are rare, and a relatively small fraction of particles on the bed surface participate in transport 
[@Ancey2010; @Roseberry2012; @Fathel2015; @Wu2019]. For example, laboratory experiments 
using a downstream light table counting device and conducted at roughly twice the threshold for particle motion involve 
the transport of less than approximately 12% of particles on the upstream bed surface [@Chartrand2017]. This result 
highlights that the flux measured across a boundary or within an area of bed surface is directly linked to the motions 
of individual particles arriving from upstream locations [@Furbish2012]. 

Because particle movements are controlled by fluid turbulence, the irregular bed surface, and collective movement effects 
[@Ancey2006; @Ancey2008; @LeeJerol2018], the connection between particle movements and the bedload 
transport rate has been difficult to formulate mathematically. The pySBeLT model provides an extensible framework within 
Python to numerically examine correlations between upstream particle entrainment rates and travel distances, with downstream 
flux. The pySBeLT model was motivated by a birth-death, immigration-emigration Markov model for bedload transport. Here, 
the movements of individual particles are represented by stochastic entrainment, motion, and deposition processes, and sediment 
flux is represented as a counting phenomenon where the number of particles in motion above the bed surface is a random 
variable [@Ancey2008]. The pySBeLT model supports ensemble simulations so that repeat numerical experiments can be conducted,
or the problem can be efficiently probed across a range of input parameter values (discussed below).

pyBeLT is run forward in time according to default or user specified parameter values in **'param.yaml'** (see the README.md for 
more details). After initialization, pyBeLT first constructs a bed of fixed particles of set_diam in both the downstream and 
cross-stream dimensions (one particle wide in the present build), and over a downstream domain length **'x_max'**. Bed surface particles
of **'set_diam'** are then randomly placed at vertices between fixed bed particles until the **'pack_density'** is met. Vertices are defined 
by a contact point between two adjacent particles. The bed of surface particles is then separated into **'num_subregions"**, and at this 
point the forward simulations are ready to commence. Simulation iterations involve three steps: (1) the number of particle entrainment 
events per **'num_subregions'** are drawn from a Poisson pmf, and this is done randomly for each numerical step up to **'n_iterations'**; 
(2) surface particles from each subregion are randomly selected for entrainment, and if there are insufficient surface particles available 
for entrainment, then all available particles are moved; (3) each entrained particle moves a distance according to a randomly sampled value 
from either the normal or lognormal distribution, and is placed at the nearest vertex between two particles that is available for placement. 
Placed particles are permitted to stack up to the **'level_limit in height'**. Travel distances of particles that exceed **'x_max'** are returned and 
queued at the upstream boundary, and are introduced back into the domain at the next numerical step according to travel distance sampling 
described above. This overall process repeats for the specified n_iterations. 

pyBeLT tracks a number of different parameters through a simulation: the vertical and horizontal positions of every particle center, 
the randomly sampled number of entrainment events, the number of particles actually entrained, the randomly sampled particle travel 
distance, the actual particle travel distance, the particle ‘age’, or the number of numerical steps since last entrainment for every 
particle, and the number of particles which cross all boundaries, i.e. sub-region and downstream at x_max. All values are written to ____ files. 
The pyBeLT model produces a time varying signal of particle flux counted at the downstream domain (as well as internal subregion domains), with a 
particle bed that changes through particle stacking and pile removal, and downstream motions of travel distance (Fig. 1). An implication 
of particle stacking within the context of the pySBeLT stochastic framework is a time varying signal of the average “particle age”, as well 
as the average “particle age range”, defined as the difference of the maximum and minimum particle ages. The model can be readily modified 
to simulate kinematics using different probability distributions, or examining particle age dynamics for deeper beds of particles available 
for transport. py_BeLT can also be extended to 2-dimensions in the cross-stream. The relatively simple parameterization of pySBeLT execution 
also makes it suitable for use as a teaching tool within advanced undergraduate and graduate courses emphasizing bedload transport.

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

@article{Ancey2006,
	title = {Statistical description of sediment transport experiments},
	volume = {74},
	url = {https://link.aps.org/doi/10.1103/PhysRevE.74.011302},
	doi = {10.1103/PhysRevE.74.011302},
	number = {1},
	journal = {Physical Review E},
	author = {Ancey, Christophe and Böhm, Tobias and Jodeau, Magali and Frey, Philippe},
	month = jul,
	year = {2006},
	note = {Publisher: American Physical Society},
	pages = {11302--11302}}
  
@article{Ancey2008,
	title = {Entrainment and motion of coarse particles in a shallow water stream down a steep slope},
	volume = {595},
	doi = {10.1017/S0022112007008774},
	journal = {Journal of Fluid Mechanics},
	author = {Ancey, C. and Davison, A. C. and Böhm, T and Jodeau, M and Frey, P},
	year = {2008},
	note = {Publisher: Cambridge University Press},
	pages = {83--114}}
  
@article{Ancey2010,
	title = {Stochastic modeling in sediment dynamics: {Exner} equation for planar bed incipient bed load transport conditions},
	volume = {115},
	url = {http://doi.wiley.com/10.1029/2009JF001260},
	doi = {10.1029/2009JF001260},
	number = {F2},
	journal = {Journal of Geophysical Research: Earth Surface},
	author = {Ancey, Christophe},
	month = jun,
	year = {2010},
	keywords = {sediment transport, Exner equation, stochastic differential equation}}
  
@article{Ancey2020,
	title = {Bedload transport: a walk between randomness and determinism. {Part} 1. {The} state of the art},
	volume = {58},
	issn = {0022-1686},
	url = {https://doi.org/10.1080/00221686.2019.1702594},
	doi = {10.1080/00221686.2019.1702594},
	number = {1},
	journal = {Journal of Hydraulic Research},
	author = {Ancey, Christophe},
	month = jan,
	year = {2020},
	note = {Publisher: Taylor \& Francis},
	pages = {1--17}}

@phdthesis{Chartrand_2017, 
	series={Electronic Theses and Dissertations (ETDs) 2008+}, 
	title={Pool-riffle dynamics in mountain streams : implications for maintenance, formation and equilibrium}, 
	url={https://open.library.ubc.ca/collections/ubctheses/24/items/1.0349138}, DOI={http://dx.doi.org/10.14288/1.0349138}, 
	school={University of British Columbia}, 
	author={Chartrand, Shawn M.}, 
	year={2017}, 
	collection={Electronic Theses and Dissertations (ETDs) 2008+}}

@phdthesis{Einstein1937,
	title = {Bedload transport as a probability problem},
	author = {Einstein, H.A.},
	year = {1937},
	note = {Pages: 53}}
  
@article{Fathel2015,
	title = {Experimental evidence of statistical ensemble behavior in bed load sediment transport},
	volume = {120},
	url = {http://doi.wiley.com/10.1002/2015JF003552},
	doi = {10.1002/2015JF003552},
	number = {11},
	journal = {Journal of Geophysical Research: Earth Surface},
	author = {Fathel, Siobhan L. and Furbish, David Jon and Schmeeckle, Mark W.},
	month = nov,
	year = {2015},
	pages = {2298--2317}}
  
@article{Furbish2012,
	title = {A probabilistic description of the bed load sediment flux: 1. {Theory}},
	volume = {117},
	url = {http://dx.doi.org/10.1029/2012JF002352},
	doi = {10.1029/2012JF002352},
	number = {F3},
	journal = {Journal of Geophysical Research: Earth Surface},
	author = {Furbish, David Jon and Haff, Peter K and Roseberry, John C and Schmeeckle, Mark W},
	month = sep,
	year = {2012},
	keywords = {bed load sediment, Fokker-Planck equation, 1825 Geomorphology: fluvial, 1862 Sediment transport, advection, diffusion},
	pages = {F03031--F03031}}
  
@article{Furb2016,
	title = {Probability distributions of bed load particle velocities, accelerations, hop distances, and travel times informed by {Jaynes}'s principle of maximum entropy},
	volume = {121},
	url = {http://doi.wiley.com/10.1002/2016JF003833},
	doi = {10.1002/2016JF003833},
	number = {7},
	journal = {Journal of Geophysical Research: Earth Surface},
	author = {Furbish, David Jon and Schmeeckle, Mark W. and Schumer, Rina and Fathel, Siobhan L.},
	month = jul,
	year = {2016},
	pages = {1373--1390}}

@article{FurbDoane2021,
	title = {Rarefied particle motions on hillslopes – {Part} 4: {Philosophy}},
	volume = {9},
	url = {https://esurf.copernicus.org/articles/9/629/2021/},
	doi = {10.5194/esurf-9-629-2021},
	number = {3},
	journal = {Earth Surface Dynamics},
	author = {Furbish, D. J. and Doane, T. H.},
	year = {2021},
	pages = {629--664}}
  
@article{LeeJerol2018,
	title = {Determining the scales of collective entrainment in collision-driven bed load},
	volume = {6},
	url = {https://www.earth-surf-dynam.net/6/1089/2018/},
	doi = {10.5194/esurf-6-1089-2018},
	number = {4},
	journal = {Earth Surface Dynamics},
	author = {Lee, D B and Jerolmack, D},
	year = {2018},
	pages = {1089--1099}}

@incollection{Parker2008, 
	address = {Reston, VA}, 
	title = {Transport of gravel and sediment mixtures}, 
	booktitle = {Sedimentation Engineering: Theory, Measurements, Modeling and Practice (ASCE Manuals and Reports on Engineering Practice No. 110)}, 
	publisher = {ASCE}, 
	author = {Parker, Gary}, 
	editor = {Garcia, MH}, 
	year = {2008}, 
	pages = {165--251}}

@article{Roseberry2012,
	title = {A probabilistic description of the bed load sediment flux: 2. {Particle} activity and motions},
	volume = {117},
	url = {http://doi.wiley.com/10.1029/2012JF002353},
	doi = {10.1029/2012JF002353},
	number = {F3},
	journal = {Journal of Geophysical Research: Earth Surface},
	author = {Roseberry, John C. and Schmeeckle, Mark W. and Furbish, David Jon},
	month = sep,
	year = {2012},
	pages = {F03032--F03032}}
  
@article{Wainwright2015,
	author = {Wainwright, John and Parsons, Anthony J. and Cooper, James R. and Gao, Peng and Gillies, John A. and Mao, Luca and Orford, Julian D. and Knight, Peter G.},
	title = {The concept of transport capacity in geomorphology},
	journal = {Reviews of Geophysics},
	volume = {53},
	number = {4},
	pages = {1155-1202},
	keywords = {sediment transport, geomorphology, turbulence, complex systems, models, management},
	doi = {https://doi.org/10.1002/2014RG000474},
	year = {2015}}

@article{Wu2020,
	title = {Generalization of {Hop} {Distance}-{Time} {Scaling} and {Particle} {Velocity} {Distributions} via a {Two}-{Regime} {Formalism} of {Bedload} {Particle} {Motions}},
	volume = {56},
	issn = {0043-1397},
	url = {https://doi.org/10.1029/2019WR025116},
	doi = {10.1029/2019WR025116},
	number = {1},
	urldate = {2022-01-13},
	journal = {Water Resources Research},
	author = {Wu, Zi and Furbish, David and Foufoula-Georgiou, Efi},
	month = jan,
	year = {2020},
	note = {Publisher: John Wiley \& Sons, Ltd},
	pages = {e2019WR025116}}
