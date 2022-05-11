## Theory Background

Sediment particle transport by **'py_SBeLT'** is simplified from natural rivers in two key ways. First, the model simulates transport along a profile 
that is one particle wide. Second, the causes of particle motions (e.g. local velocities or shear stresses) are represented in the model through the 
resulting kinematics--namely the number of entrainment events and the particle hop or travel distance. These two quantities alone are the main components 
that drive transport and the resulting behavior that emerges from model simulations with **'py_SBeLT'**. The remainder of the code is primarily keeping 
track of particle motions, the changing bed conditions and locations along the bed that are available for deposition. 

### Entrainment

With this general set-up, transport in **'py_SBeLT'** is treated as a stochastic process. The number of entrainment events and the particle travel 
distance are random variables that are indexed to time and space. Entrainment events are sampled from the Poisson probability mass function (pmf) (Fig. 
1). The Poisson pmf expresses the probability that a specific number of events will occur within a time interval according to the rate constant &#955;. 
This means particle entrainment is treated as independent events between the **'num_subregions'** (see README.md and paper.md) and between each 
iteration. Therefore, the Poisson pmf is fixed for any given simulation (see DEFAULT_PARAMS.md). However, for each iteration and subregion a new value is 
randomly sampled from the Poisson pmf to define the entrainment events. 

|![Image](/figures/poisson.png)
|:--:| 
| *Figure 1. Example Poisson pmf used to set the number of entrainment events in **`py_SBeLT`**.* |

Use of the Poisson pmf is linked to sediment transport theory and supporting experimental results which show **'rarefied'** transport (Furbish et al., 
2016) can be characterized as a Poisson process when entrainment includes effects related to fluid phenomena under steady state transport conditions
(Ancey et al., 2008). Collective related entrainment effects are not represented (Ancey et al., 2008; Lee and Jerolmack, 2018). We have tested 
2018) **'py_SBeLT'** for a range of rate constant values (see readme.md). Results from our testing reveal that the value specified for &#955;<sub>1</sub> 
2019) along with the **'num_subregions'** controls the intensity, or magnitude of transport. This makes intuitive sense because the local transport 
2020) intensity is directly linked to the number of sediment particles entrained from the bed surface.

### Travel Distance

Physical experiments show that particle travel distances under **'rarefied'** transport conditions are commonly skewed to longer lengths with a 
well defined mode >0 that is on the order of 1 or more particle diameters (Lajueness et al., 2010; Fathel et al., 2015). Fathel et al. (2015, page 2309) 
describe this tendancy in the following manner:

> That is, particles initially have a high likelihood of disentrainment but then experience a decreasing spatial disentrainment rate...with increasing 
> hop distance L<sub>x</sub>.

**'py_SBeLT'** notably does not consider disentrainment rates because entrainment and deposition occur within the same numerical iteration. However, the 
general idea expressed by Fathel et al. (2015) is that following entrainment, sediment particles are most likely to deposit nearest to the point of 
entrainment, with decreasing liklihood at points increasingly distant from the location of entrainment. For **'py_SBeLT'** this is entirely controled by 
the randomly sampled travel distance and the availability of deposition locations in proximity to the mapped travel distance. 

The particle bed in the model is characterized by a discrete number of available deposition locations because these locations are defined by vertices 
between adjacent particles (subject to the level limit--see DEFAULT_PARAMS.md). As a result, the total number of deposition locations scales 
with the total number of particles which define the domain length. This may appear problematic because, in essence, this model setup may bias simulated 
particle motions by limiting travel distances to a finite number of deposition locations. However, sediment particle disentrainment in rivers is 
influenced by the position of particles on the bed surface, and previous publications commonly express travel distance as a probability mass function, 
or a discrete density function (Lajueness et al., 2010; Fathel et al., 2015). 

Model testing revealed two important outcomes with respect to parameterizing the particle travel distance, noting that **'py_SBeLT'** samples distances 
from probability density functions (pdf) because it simplifies the process of finding available deposition sites. Travel distance probability density 
functions with modes increasingly close to zero, or scaling as a few particle diameter equivalents in length or less, leads to two challenges. First, the 
development of isolated particle piles, and second, realized particle travel distances that increasingly diverge from the underlying randmoly sampled 
distance. For example, if the travel distance mode is 1 or 2 particle diameter equivalents in length, the model searches for deposition locations which 
are relatively close to the point of entrainment. When this occurs at multiple locations on the bed, over a relatively short period of time particle
piles emerge--i.e. particles stack and form triangular piles (see gif at github home page). 

Because **'py_SBeLT'** simulates rarefied transport conditions, these piles tend to initially occupy isolated posiitons in space. However, as the 
particle piles grow in height up to a vertical height limit (see readme.md and paper.md), locally available deposition locations are increasingly distant 
from entrainment locations. Under this scenario, future entrained particles travel to more distant locations regardless of their sampled travel distance 
because the liklihood is high that deposition sites relatively close to the point of entrainment will not be available. The model is therefore forced to 
search farther beyond the sampled travel distance for an available depoistion location. This leads to enlongated growth of the initially isoloated piles, 
and anomalous transport behavior. 

We overcame these two challenges by using probability density functions which provide for modes displaced from zero (generally > 1 or more particle 
diameter equivalents in length), and for which the probability of sampling relatively small values Pr(X<=x) vanishes as x &#8594; 0. At present, 
**'py_SBeLT'** uses the normal or lognormal probability density functions to specify sediment particle travel distances (Fig. 2). The normal distribution 
notably does not satisfy observations of travel distance distributions that are skewed to long lengths, and it can also yield negative travel distances.
With respect to the first issue, we have incorporated the normal distribution within the model build as a reference condition for generating
quasi-random sets of travel distances, consistent with numerous statistical sampling approaches. We address the second issue by limiting possible travel
distance values to >= 0 when the normal distribution is selected.

|![Image](/figures/lognormal.png)
|:--:| 
| *Figure 2. Example Lognormal pdf used to set the particle travel distance in **`py_SBeLT`**.* |

The selected probability density function is specified by the user or default settings can be used (see README.md, paper.md and DEFAULT_PARAMS.md).
The readme.md provides the distribution function parameter value ranges tested to date. Additional probability density functions can also be used to 
sample particle travel distances. For example, the gamma distribution can provide modes displaced from zero, and with distribution shapes that are 
skewed to longer lengths. Physical experiments that in part motivated development of **`py_SBeLT`** report that particle hop distances are 
well described by the Weibull distribution (Fathel et al., 2015)(Fig. 3). Note however that **'py_SBeLT'** would likely produce isolated particle piles 
if the shape parameter *k* of the Weibull distribution were constrained to values <1 (Fathel et al., 2015). This limitation reflects the present model
build which simulates rarefied transport along a downstream profile of particles with available deposition locations constrained by particle diameter. We 
expect that a model build expanded to cross-stream scales of many particle diameters in length might be more capable of simulating the experimental 
conditions reported by Fathel et al. (2015). An additional probability density function that could be incorporated for travel distance sampling is the 
Pareto distribution, which is subject to a minimum value for the model random variable. In the present case the minimum value could be set to ~> 1 or 
more particle diameter equivalents in length in order to avoid the particle piling issue described above. The Pareto distribution, like the lognormal 
distribution, also has the advantage of being heavy-tailed where the variance becomes undefined for certain shape parameters. This provides an 
additional opportunity to illustrate the consequences of heavy-tailed distributions and the influence on the time-series of sediment flux. These 
additional and other pdfs can be easily incorporated into future model extensions using the **'scipy.stats'** library of functions.

|![Image](/figures/weibull.png)
|:--:| 
| *Figure 3. Example Weibull pdf that could be used to set the number of entrainment events in **`py_SBeLT`**.* |

## References

Ancey, C., Davison, A. C., Böhm, T., Jodeau, M., & Frey, P. (2008). Entrainment and motion of coarse particles in a shallow water stream down a steep 
slope. Journal of Fluid Mechanics, 595, 83–114. https://doi.org/10.1017/S0022112007008774.

Fathel, S. L., Furbish, D. J., & Schmeeckle, M. W. (2015). Experimental evidence of statistical ensemble behavior in bed load sediment transport. Journal 
of Geophysical Research: Earth Surface, 120(11), 2298–2317. https://doi.org/10.1002/2015JF003552.

Furbish, D. J., Schmeeckle, M. W., Schumer, R., & Fathel, S. L. (2016). Probability distributions of bed load particle velocities, accelerations, hop 
distances, and travel times informed by Jaynes’s principle of maximum entropy. Journal of Geophysical Research: Earth Surface, 121(7), 1373–1390. 
https://doi.org/10.1002/2016JF003833.

Lajeunesse, E., Malverti, L., & Charru, F. (2010). Bed load transport in turbulent flow at the grain scale: Experiments and modeling. Journal of 
Geophysical Research: Earth Surface, 115(F4). https://doi.org/10.1029/2009JF001628.

Lee, D. B., & Jerolmack, D. (2018). Determining the scales of collective entrainment in collision-driven bed load. Earth Surface Dynamics, 6(4), 
1089–1099. https://doi.org/10.5194/esurf-6-1089-2018.
