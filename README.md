# Py_SBeLT

## Requirements

Installing Python requirements:

```bash
 pip3 install -r requirements.txt
```

## Nomeculture

- **Bed Particle**: A particle embedded in the stream bed. Does not move.
- **Model Particle**: A particle subject to entrainment events. Moves.
- **Event Particle**: A Model Particle which has been selected for entrainment.
- **Desired Hop**: The hop distance an Event Particle is assigned at the beginning of an entrainment event
- **Stream**: Consists of the Bed and Model particles.

## Parameters

The model reads in parameters from a yaml file. The default parameters file for the model is `param.yaml`, located in the **`model/`** directory:

```
project
│   README.md
│   requirements.txt
└───model
│   │   param.yaml
│   │   run.py
│   │   ...
```

Comments indicate what values have been tested for each parameter and within what range. Users are welcome to enter parameters outside these bounds to see the results, however data type violations and specific values will be caught during the validation process and cause the model not to run (e.g set_diam: 0).

## Running the Model

First, confirm you have a **`plots/`** directory located in the **`BeRCM/`** project:

```
BeRCM
│   README.md
│   requirements.txt
└───model
│   │   ...
└───plots <--- Make sure this exists!
│   │   ...
```

### Running from Command Line

1. Set current directory to **`BeRCM/model`**:

```bash
 cd /path/to/BeRCM/model
```

2. 
    a. Run a _single_ instance of the model:

    ```bash
    python3 run.py PARAM_FILE
    ```

    Where `PARAM_FILE` is the path to the desired parameters file. To use the default parameters file set `PARAM_FILE=`**`param.yaml`**. 

    b. Run _multiple_ instances of the model:

    ```bash
    python3 parallel.py NUM_PROCESSES PARAM_FILE...
    ```

    Note that multiple parameter files can be passed to the **`parallel.py`** script. If more than one parameter file is passed then the number of files passed must be equal to the number of processes requested.

### Running in Spyder (THIS SECTION IS WIP)

<!-- 1. Open **`run.py`** and **`parameters.py`** in Spyder

2. Navigate to **`run.py`** and execute using the `'Run file (F5)'` button

If you wish to view the plots in the plot pane instead of saving, set the last call 
to `stream` and `flux_info` at the end of **`run.py`** (58-60) to:

```python
    ml.plot_stream(iteration, bed_particles, model_particles, pm.x_max, 10, 
                   available_vertices, to_file=False)
```

and:

```python
    plot.flux_info(particle_flux_list, to_file=False)
```

Ensuring `to_file` is set to `False` for one or both calls. Then follow steps 1 and 2 above. -->

## Results

All relevant information from the run will be stored in a Python [shelf](https://docs.python.org/3/library/shelve.html) in the **`plots/`** directory:

```
RCM_BedloadTransport
│   README.md
│   requirements.txt
└───model
│   │   ...
└───plots <--- here!
│   │   ...
```

Each run will have a unique shelf filename based on the following syntax: `run-info-YYYYMM-DDHH-MMSS-UUID`, where `YYYYMM-DDHH-MMSS` represents the time the model was started and `UUID` is a unique identifier assigned by the script. For example, a run started on December 24th, 2030 at 12:01:42pm and assigned the UUID 1a4c536f-4bf3-4534-becf-90f207caf9e will use the filename: 
    
- **`run-info-203012-2412-0142-1a4c536f-4bf3-4534-becf-90f207caf9e`**

Each shelf file will contain: the parameters used, the bed particle array, the model particle array, available vertices, and event particles for each iteration, the particle flux list, and the particle age and age range lists.

Each of these data values can be retrieved by opening the shelf file and using the appropriate keys as shown below:

```{python3}
with shelve.open(beRCM_shelf, 'r') as shelf:
    shelf["param"]      # parameters
    shelf["bed"]        # bed particles
    shelf["12"][0]      # model particles at end of iteration 12
    shelf["12"][1]      # available vertices for iteration 12
    shelf["12"][2]      # events particles for iteration 12
    shelf["flux"]       # flux list
    shelf["avg_age"]    # avg age list
    shelf["age_range"]  # age range list
```

## Plotting
### Plotting from Shelf
The project currently contains logic for plotting a visual of the streambed, the particle flux distribution, and particle age-related information.

Navigate to the **`plots/`** directory. Choose a desired run-info shelf file to plot, create or choose a location for the plots to save, and decide what range of iterations to plot for the streambed. Then run the following command:

```{bash}
python3 plot_maker.py SHELF_NAME SAVE_LOCATION MIN_ITER MAX_ITER
```

### Creating Gif of Streambed

A GIF can me made from streambed plots using the **`gif_maker.py`** script. In the script alter the `in_dir` and `out_dir` variables to point to the location of the strambed plots and where the GIF should be saved resepectively. Similarly, change the `start` and `stop` variables to reflect which range of plots to use in the GIF compilation.

For example:

```{python3}
in_dir = '../plots/test/'
out_dir = '../plots/test/'

start = 0
stop = 99
```

Will create a GIF using plots of iterations 0-99 and will store the GIF in the same location of the .png plots.

<!-- Embed a gif here as an example/motivation -->
