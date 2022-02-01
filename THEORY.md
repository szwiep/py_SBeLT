## Theory

Sediment particle transport in **'py_SBeLT'** is simplified from natural rivers in two key ways. Transport is simulated along a profile that is one particle wide,
and the causes of particle motions are represented through the resulting kinematic descriptions--namely the number of entrainment events and the particle hop or
travel distance. These two quatities alone are the engine that drives transport and the resulting behavior that emerges from model simulations. Most of the code
is simply keeping track of particle motions and locations along the bed that are available for deposition. 

With this general set-up, transport in **'py_SBeLT'** is treated as a stochastic process, with two random variables that are indexed to time and space: the number
of entrainment events and the particle travel distance. Entrainment events are sampled from the Poisson probability mass function (pmf). The Poisson pmf expresses
the probability that a specific number of events will occur within a time interval according to the rate constant &#955;. Practically, this means 
entrainment is treated as independent events between **'num_subregions'** (see readme.md and paper.md) and between each iteration. For example, 
$\lambda$ is specified within the **'parameter.yaml'** file as &#955~1~. Therefore, the Poisson pmf is fixed for any given 
simulation. However, for each iteration and subregion a new value is randomly sampled from the Poisson pmf to define the entrainment events. 

Use of the Poisson pmf is linked to sediment transport theory and associated physical experiments which indicate **'rarefied'** transport [@Furb:2016] is 
characterized as a Poisson process when entrainment includes effects related to fluid phenomena under steady state transport conditions [@Ancey:2008]. Collective 
related entrainment effects are not represented [@Ancey:2008; @LeeJerol:2018]. We have tested **'py_SBeLT'** for a range of rate parameter values (see readme.md). 
Results from this testing reveals that the value specified for $\lambda_1$ along with the **'num_subregions'** controls the intensity or magnitude 
of transport.

Physical experiments have generally shown that particle travel distances under **'rarefied'** transport conditions is commonly skewed to longer lengths with a 
well defined mode $>$0, and on the order of 1 or more particle diameters (Lajueness et al., 2010; Fathel et al., 2015). Fathel et al. (2015) describe this 
tendancy in the following manner:

> That is, particles initially have a high likelihood of disentrainment but then experience a decreasing spatial disentrainment rate...with increasing hop distance L~x~.

**'py_SBeLT'** does not consider disentrainment rates because entrainment and deposition occur within the same numerical iteration. However, the general idea 
expressed by Fathel et al. (2015) is that following entrainment, sediment particles are most likely to deposit nearest to the point of entrainment, with 
decreasing liklihood at points increasingly distant from the location of entrainment. For **'py_SBeLT'** this is entirely controled by the randomly sampled travel 
distance and the availability of deposition locations in proximity to the mapped travel distance. The particle bed in the model is characterized by a discrete 
distribution of available deposition locations because these locations are defined by vertices between adjacent particles (subject to the level limit specified in 
**'parameter.yaml'**). 

Due in part to the simplicity of the model build, testing revealed two important outcomes with respect to parameterizing the particel travel distance. Travel
distance distributions with modes close to zero lead to: (1) the development of isolated particle piles, and (2) realized particle travel distances that diverge 
from the underlying randmoly sampled distance. For example, if isolated particle piles occur and are limited by the height limit (see readme.md and paper.md), 
locally available deposition locations are increasingly distant from the entrainment location. Therefore, for travel distance modes close to zero, the closest 
available deposition location has a high liklihood of diverging from the nearest location associated with the randomly sampled travel distance.

We overcame this challenge in the most reasonable manner possible by using probability distribution functions which provide for modes displaced from zero, and for 
which the probability of sampling relatively small values vanishes $\rightarrow$0 

