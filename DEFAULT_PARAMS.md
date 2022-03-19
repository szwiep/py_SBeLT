## Default Model Parameter Values and Settings

Here we briefly describe the default model parameter values and settings for **'py_SBeLT'**. Readers are directed to **'README.md'** for information
concerning the ranges of parameter values that have been testing using **'py_SBeLT'**.

### Iterations

**Default Value = 1000**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Bed_length

**Default Value = 100**

We use a default value of 100 for *Bed_length*. This length of bed was found to provide *sufficient* space to observe transport behaviors that are common
to natural rivers, such as variable sediment particle residence times in the bed, spatially variable rates of transport, etc. The model was tested up to
a *Bed_length* of 1000 *millimeters* are performed as expected. It is important to note however that the combination of *Bed_length* and *Num_subregions*
in conjunction with *Poiss_lambda* control the intensity of transport simulated by **'py_SBeLT'**. Therefore, if users change the parameter *Bed_length*,
it may be useful to also consider changing the default values for *Num_subregions* and *Poiss_lambda*. 

### Particle_diam

**Default Value = 0.5**

We use a default value of 0.5 for *Particle_diam*, with assumed units of length of *millimeters*. This default value provides a reasonable scale--
assuming a default *bed_length* of 100 *millimeters*--to view simulation results when plotted over the entire domain. This default value is also
consistent with the grain sizes present on the bed surface of many low land rivers. 

### Particle_pack_den

**Default Value = 0.78**

We use a default value of 0.78 for *Particle_pack_den*. This provides a surface coverage of 78% over the *Bed_length*. More or fewer particles on the 
surface can affect derived simulated behaviors such as how long particles sit in the mobile part of bed surface before mobilization occurs, or the 
availability of deposition locations in relation to the relative particle level, or stacking height.

### Num_subregions

**Default Value = 4**

We use a default value of 4 for *Num_subregions*, and in conjuction with the default *Bed_length* of 100 *millimeters*. If the user holds *Bed_length* 
and *Poiss_lambda* to the default values, increasing the *Num_subregions* will result in increased rates of particle transport up to the limit imposed by 
*Particle_pack_den*, and decreasing the *Num_subregions* will result in less transport.

### Level_limit

**Default Value = 3**

We use a default value of 3 for *Level_limit*. This means that particles can stack 3 high in any region of the bed subject to particle availability, 
selected hop distances, etc. Using this default value provides for development of bed topography. As a result, model particles develop a unique history 
of age, or the duration between mobilization events. Using smaller values for *Level_limit* will depress bed topography, and truncate the distribution of 
particle ages.

### Poiss_lambda

**Default Value = 5**

We use a default value of 5 for *Poiss_lambda* to provide a relatively high rate of entrainment provided the default settings for *Bed_length*, 
*Particle_pack_den* and *Num_subregions*. Notably the value used for *Poiss_lambda* does not affect model behavior for the range of values tested for 
**'py_SBeLT'** (see README.md), but it does influence how many iterations may be needed to achieve a relatively high number of total entrainment events, 
or particle passing events beyond internal and the downstream model boundary. These outcomes will influence interpretations of model results.   

### Gauss

**Default Value = 'False'**

We use a default value of *False* for *Guass*. This means particle travel distances are sampled from a *Lognormal* continuous distribution, as opposed to 
a *Normal* distribution. We use the setting of *False* because this permits sampling from a distribution that can be skewed to longer distances, which is 
more consistent with experimental observations (see THEORY.md). Testing reveals however that use of a *Normal* distribution can yield behavior generally 
consistent with experimental results provided that values are truncated at zero and the mean, or central tendency, is positively offset from zero.

### Gauss_mu

**Default Value = 1**

We use a default value of 1 for *Gauss_mu*. Depending on values selected for *Gauss_sigma*, a value of 1 for *Gauss_mu* provides a range of particle 
travel distances which scale from ~1 to 20 particle diameters (for the default *Particle_diam* setting), or more (accounting for the log transformation 
of *X*). More specifically, the value selected for *Gauss_mu* sets the median value for the distribution of all possible travel distances; in the default 
case this means that half of all possible sampled travel distances will scale as more than 5 particle diameters (accounting for the log transformation 
of *X*). Testing with **'py_SBeLT'** reveals that this default setting generally avoids the development of isolated particle piles and model behavior 
which differs from experimental observations (see THEORY.md). 

### Gauss_sigma

**Default Value = 0.25**

We use a default value of 0.25 for *Gauss_sigma*. For the default travel distance distribution setting (see *Gauss*), relatively large values for 
*Gauss_sigma* provide an increasingly skewed distribtuion to larger travel distances, and small values provide a distribution which tends to a 
*Gaussian*.

### Data_save_interval

**Default Value = 1**

We use a default value of 1 for *Data_save_interval*. This default value will save model output every iteration.

### Height_dependent_entr

**Default Value = 'False'**

We use a default value of *False* for *Height_dependent_entr*. This setting means that particles sitting at the height limit set by *Level_limit* are not 
automatically entrained during the next iteration. We use a default value of *False* because a setting of *True* removes the role of entrainment events 
variability in simulated behavior of transport. We provide the option to set this parameter to *True* to permit exploring how automatic entrainment may 
influence particle age distributions, or bed topography.

### Out_path

**Default Value = '.'**

We use a default setting of writing output to the active directory.

### Out_name

**Default Value = 'sbelt-out'**

We use a default setting of writing output to *sbelt-out*.
Simulations are run using *Iterations = 10$^6$*.
