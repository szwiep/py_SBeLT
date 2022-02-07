## Theory

Sediment particle transport in **'py_SBeLT'** is simplified from natural rivers in two key ways. Transport is simulated along a profile that is one particle wide,
and the causes of particle motions are represented through the resulting kinematics--namely the number of entrainment events and the particle hop or
travel distance. These two quantities alone are the engine that drives transport and the resulting behavior that emerges from model simulations with
**'py_SBeLT'**. The remainder of the code is primarily keeping track of particle motions, the changing bed conditions and locations along the bed that are
available for deposition. 

With this general set-up, transport in **'py_SBeLT'** is treated as a stochastic process. There are two random variables that are indexed to time and space: the
number of entrainment events and the particle travel distance. Entrainment events are sampled from the Poisson probability mass function (pmf). The Poisson pmf
expresses the probability that a specific number of events will occur within a time interval according to the rate constant &#955;. Practically, this means 
particle entrainment is treated as independent events between the of **'num_subregions'** (see readme.md and paper.md) and between each iteration. For example, 
&#955; is specified within the **'parameter.yaml'** file as &#955;<sub>1</sub>. Therefore, the Poisson pmf is fixed for any given simulation. However, for each
iteration and subregion a new value is randomly sampled from the Poisson pmf to define the entrainment events. 

Use of the Poisson pmf is linked to sediment transport theory and associated physical experiments which indicate **'rarefied'** transport (Furbish et al., 2016)
is characterized as a Poisson process when entrainment includes effects related to fluid phenomena under steady state transport conditions (Ancey et al., 2008).
Collective related entrainment effects are not represented (Ancey et al., 2008; Lee and Jerolmack, 2018). We have tested **'py_SBeLT'** for a range of rate 
parameter values (see readme.md). Results from this testing reveals that the value specified for &#955;<sub>1</sub> along with the **'num_subregions'** controls
the intensity or magnitude of transport.

Physical experiments have generally shown that particle travel distances under **'rarefied'** transport conditions are commonly skewed to longer lengths with a 
well defined mode >0 that is on the order of 1 or more particle diameters (Lajueness et al., 2010; Fathel et al., 2015). Fathel et al. (2015) describe this 
tendancy in the following manner:

> That is, particles initially have a high likelihood of disentrainment but then experience a decreasing spatial disentrainment rate...with increasing hop distance L<sub>x</sub>.

**'py_SBeLT'** does not consider disentrainment rates because entrainment and deposition occur within the same numerical iteration. However, the general idea 
expressed by Fathel et al. (2015) is that following entrainment, sediment particles are most likely to deposit nearest to the point of entrainment, with 
decreasing liklihood at points increasingly distant from the location of entrainment. For **'py_SBeLT'** this is entirely controled by the randomly sampled travel 
distance and the availability of deposition locations in proximity to the mapped travel distance. The particle bed in the model is characterized by a discrete 
distribution of available deposition locations because these locations are defined by vertices between adjacent particles (subject to the level limit specified in 
**'parameter.yaml'**). 

Testing revealed two important outcomes with respect to parameterizing the particle travel distance. Travel distance distributions with modes close to zero lead
to two challenges. First, the development of isolated particle piles, and second, realized particle travel distances that diverge from the underlying randmoly
sampled distance. For example, if isolated particle piles occur and their vertical height is controlled by the height limit (see readme.md and paper.md), locally
available deposition locations are increasingly distant from the entrainment location. This is because most sediment particles are sampled to travel a relativley
short distance, and therefore the model looks for an available deposition location close to the point of entrainment. However, as more iterations occur the 
liklihood is high that deposition sites relatively close to the point of entrainment will not be available. The model is therefore forced to search farther beyond 
the sampled travel distance for an available depoistion location. This leads to the development of enlongated and isoloated piles, and anomalous transport
behavior. 

We overcame thses two challenges by using probability distribution functions which provide for modes displaced from zero (generally > 1 particle diameter equivalent), and for which the probability of sampling relatively small values Pr(X<=x) vanishes as x &#8594; 0. 

## References

Ancey, C., Davison, A. C., Böhm, T., Jodeau, M., & Frey, P. (2008). Entrainment and motion of coarse particles in a shallow water stream down a steep slope. Journal of Fluid Mechanics, 595, 83–114. https://doi.org/10.1017/S0022112007008774

Fathel, S. L., Furbish, D. J., & Schmeeckle, M. W. (2015). Experimental evidence of statistical ensemble behavior in bed load sediment transport. Journal of Geophysical Research: Earth Surface, 120(11), 2298–2317. https://doi.org/10.1002/2015JF003552

Furbish, D. J., Schmeeckle, M. W., Schumer, R., & Fathel, S. L. (2016). Probability distributions of bed load particle velocities, accelerations, hop distances, and travel times informed by Jaynes’s principle of maximum entropy. Journal of Geophysical Research: Earth Surface, 121(7), 1373–1390. https://doi.org/10.1002/2016JF003833

Lajeunesse, E., Malverti, L., & Charru, F. (2010). Bed load transport in turbulent flow at the grain scale: Experiments and modeling. Journal of Geophysical Research: Earth Surface, 115(F4). https://doi.org/10.1029/2009JF001628

Lee, D. B., & Jerolmack, D. (2018). Determining the scales of collective entrainment in collision-driven bed load. Earth Surface Dynamics, 6(4), 1089–1099. https://doi.org/10.5194/esurf-6-1089-2018



