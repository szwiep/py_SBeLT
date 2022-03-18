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
to natural rivers, such as variable sediment particle residence times in the bed, spatially variable rates of transport, etc. The model was tested up to a 
*Bed_length* of 1000 *millimeters* are performed as expected. It is important to note however that the combination of *Bed_length* and *Num_subregions* in 
conjunction with *Poiss_lambda* control the intensity of transport simulated by **'py_SBeLT'**. Therefore, if users change the parameter *Bed_length*, it 
may be useful to also consider changing the default values for *Num_subregions* and *Poiss_lambda*. 

### Particle_diam

**Default Value = 0.5**

We use a default value of 0.5 for *Particle_diam*, with assumed units of length of *millimeters*. This default value provides a reasonable scale--assuming
a default *bed_length* of 100 *millimeters*--to view simulation results when plotted over the entire domain. This default value is also consistent with the
grain sizes present on the bed surface of many low land rivers. 

### Particle_pack_den

**Default Value = 0.78**

We use a default value of 0.78 for *Particle_pack_den*. This provides a surface coverage of 78% over the *Bed_length*. More or fewer particles on the surface
can affect derived simulated behaviors such as how long particles sit in the mobile part of bed surface before mobilization occurs, or the availability of
deposition locations in relation to the relative particle level, or stacking height.

### Num_subregions

**Default Value = 4**

We use a default value of 4 for *Num_subregions*, and in conjuction with the default *Bed_length* of 100 *millimeters*. If the user holds *Bed_length* and 
*Poiss_lambda* to the default values, increasing the *Num_subregions* will result in increased rates of particle transport up to the limit imposed by *Particle_pack_den*,
and decreasing the *Num_subregions* will result in less transport.

### Level_limit

**Default Value = 3**

We use a default value of 3 for *Level_limit*. This means that particles can stack 3 high in any region of the bed subject to particle availability, selected
hop distances, etc. Using this default value provides for development of bed topography. As a result, model particles develop a unique history of age, or 
the duration between mobilization events. Using smaller values for *Level_limit* will depress bed topography, and truncate the distribution of particle ages.

### Poiss_lambda

**Default Value = 5**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Gauss

**Default Value = 'False'**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Gauss_mu

**Default Value = 1**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Gauss_sigma

**Default Value = 0.25**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Data_save_interval

**Default Value = 1**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Height_dependent_entr

**Default Value = 'False'**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Out_path

**Default Value = '.'**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.

### Out_name

**Default Value = 'sbelt-out'**

We use a default value of 1000 for *Iterations* to allow users to explore **'py_SBeLT'** in an efficient way without having to wait to long for results.
Simulations are run using *Iterations = 10$^6$*.
