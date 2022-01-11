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

Granular sediment of various sizes moves downstream along river beds when water 
flow is capable of entraining particles from the bed surface. This process is known 
as bedload sediment transport because the particles travel close to the boundary. 
It is common to treat the transport process as a predictive problem, in which the 
mean transport rate of particles at a stationary observation point is a function of 
local water flow conditions (Parker, 2008; Wainwright, 2014; Ancey, 2020). Testing 
transport functions against data from natural rivers reveals generally inaccurate 
results (Wilcock, 2001; Ancey, 2020a). Calibration of transport functions to local 
river conditions can provide improved results (Wilcock, 2001). However, treating 
bedload sediment transport as a deterministic phenomenon neglects the stochastic 
nature of transport (Einstein, 1937; Ancey et al., 2006; Furbish, 2012; 
Wainwright 2014; Fathel et al., 2015; Ancey 2020a, Pierce and Hassan, 2020). 
Here, we present an open-source Python software, pyBeLT, which simulates the 
kinematics of rarefied particle transport (low rates) as a stochastic process along 
a riverbed profile.

# Statement of need

Research at the intersection of geomorphology, geophysics and hydraulics is increasingly 
focused on building a theoretical foundation for the treatment of bedload transport as a 
stochastic phenomenon–building from ideas introduced by Einstein (1937, 1950) 
(Ancey et al., 2006; Ancey et al., 2008; Ancey, 2010; Furbish and Haff, 2010; Furbish, 2012; 
Furbish and Schmeeckle, 2013; Fathel et al., 2015; Pierce and Hassan, 2020; Ancey, 2020a, 2020b; 
Furbish and Doane, 2021). Associated theories are commonly tested against data derived from 
laboratory experiments (e.g. Ancey, 2010; Roseberry et al., 2012; Fathel et al., 2015; Wu et al., 2019), 
in which transport magnitude is low to moderate, reflecting rarefied conditions-i.e. noncontinuum 
(Furbish et al., 2016). Under rarefied transport interactions between two or more particles while in 
motion are negligible (Furbish and Doane, 2021), and a relatively small fraction of particles on the 
bed surface participate in transport. For example, laboratory experiments conducted at roughly a 
factor 1.5-2 above the threshold flow conditions for particle entrainment included the transport of 
approximately <1-12%, on average, of all particles on a unit area of the bed surface measuring the 
local channel width in both dimensions (Chartrand, 2017).

Laboratory experiments and associated data used to test stochastic transport theories involves 
tracking individual particle motions within relatively small sampling windows commonly measuring 
order tens or more particle diameters in downstream length (Ancey et al., 2006; Roseberry et al., 2012;
Fathel et al., 2015). This introduces challenges related to bed load transport information censorship 
(Ballio et al., 2019), as well as the software required to track particle motions, and the relatively 
long processing times required to extract kinematic information. Here, we report software that partially 
overcomes these challenges in the most general way we can imagine, within a Python framework that supports 
ensemble numerical experiments. The software pyBeLT was motivated by a birth-death, immigration-emigration 
Markov model for bed load transport. Hence, bed load transport is treated as a counting phenomenon where 
the number of particles in motion N within a control area A above the bed surface is a random variable 
(Ancey et al., 2008). Below we review the basic steps followed by pyBeLT to simulate bed load transport 
in rivers under rarefied conditions.

pyBeLT is run forward in time according to default or user specified parameter values in param.yaml: 
the initial surface packing density of placed particles pack_density, the length of the domain x_max, 
the particle diameter set_diam, the number of bed surface sub-regions num_subregions, the number of 
vertical particle stacking levels level_limit (positive integer =>1), the number of numerical steps n_iterations,
the Poisson pmf rate constant lambda which sets the number of entrainment events En per numerical step and 
for each sub-region lambda_1, a boolean to specify the probability distribution which sets the entrained 
particle travel distance (true-normal distribution, or false-lognormal distribution), expected value of the 
travel distance distribution mu, the standard deviation of the travel distance distribution sigma, the interval 
between model output saves snapshot_interval and a boolean setting the state of height dependent entrainment of 
particles (true-particles occurring in the highest positions are entrained at the next simulation step, or 
false-particles occurring in the highest positions entrain only if randomly selected according to lambda_1.

After initialization, pyBeLT first constructs a bed of fixed particles of set_diam in both the downstream and 
cross-stream dimensions, and over a downstream domain length x_max. Bed surface particles of set_diam are then 
randomly placed at vertices between fixed bed particles until the pack_density is met. Vertices are defined by 
a contact point between two adjacent particles. The bed of surface particles is then separated into num_subregions, 
and at this point the model is set-up to perform simulations. Simulation iterations involve three steps: (1) the 
number of particle entrainment events En per num_subregions are drawn from a Poisson pmf, and this is done randomly 
for each numerical step up to n_iterations; (2) En surface particles from each subregion are randomly selected for 
entrainment, and if there are insufficient surface particles available for entrainment, then all available particles 
are moved; (3) each entrained particle moves a distance according to a randomly sampled value from either the normal 
or lognormal distribution, and is placed at the nearest vertex between two particles that is available for placement. 
Placed particles are permitted to stack up to the level_limit in height. Travel distances of particles that exceed x_max 
are returned and queued at the upstream model domain, and are introduced back into the model at the next numerical step 
according to travel distance sampling described above. This overall process repeats for the specified n_iterations. 

pyBeLT tracks a number of different parameters through a simulation: the vertical and horizontal positions of every 
particle center, the randomly sampled number of entrainment events En, the number of particles actually entrained, 
the randomly sampled particle travel distance, the actual particle travel distance, the particle ‘age’, or the number of 
numerical steps since last entrainment for every particle, and the number of particles which cross all model boundaries–
sub-region and downstream at x_max. All values are written to ____ files. Ensemble simulations using the same param.yaml 
inputs, or multiple different param.yaml inputs is possible (see README.md). pyBeLT produces a time varying signal of 
particle flux counted at the downstream domain (as well as internal subregion domains), with a particle bed that changes 
through particle stacking and pile removal, and downstream motions of travel distance \autoref{fig:fig1}. An implication 
of particle stacking within the context of a stochastic model framework is a time varying signal of the average “particle age”, 
as well as the average “particle age range”, defined as the difference of the maximum and minimum particle ages. The software 
can be readily modified to simulate kinematics using different probability distributions, or examining particle age dynamics 
for deeper beds of particles available for transport. py_BeLT can also be extended to 2-dimensions in the cross-stream. 
The relatively simple parameterization of model execution also makes it suitable for use as a teaching tool within advanced 
undergraduate and graduate courses emphasizing bedload transport.

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

<p>
![Figure 1](/figures/Figure 1.pdf)
</p>

# Acknowledgements

S. Zwiep was funded in part through an Undergraduate Student Research Award from the 
National Science and Engineering Research Council of Canada (NSERC). S.M. Chartrand was 
funded through a Postdoctoral Fellowship awarded by NSERC, and through internal research 
funding provided by Simon Fraser University. The model was inspired by discussions with 
David Jon Furbish, who also provided useful input and critical feedback at various stages 
of model development and testing. Kevin Pierce also provided helpful feedback during model 
development. 


# References
