# Py_BeRCM

### Requirements

Installing Python requirements:

```bash
 pip3 install -r requirements.txt
```

### Nomeculture 

- **Bed Particle**: A particle embedded in the stream bed. Does not move.
- **Model Particle**: A particle subject to entrainment events. Moves.
- **Event Particle**: A Model Particle which has been selected for entrainment.
- **Desired Hop**: The hop distance an Event Particle is assigned at the beginning of an entrainment event
- **Stream**: Consists of the Bed and Model particles.

### Parameters

The model comes with a default parameters file, `parameters.yaml`, which users can edit. The file is located in the **`model`** directory:

```
project
│   README.md
│   requirements.txt
└───model
│   │   parameters.py
│   │   run.py
│   │   ...
```

Comments indicate what values have been tested for each parameter and within what range. Users are welcome to enter parameters outside these bounds to see the results, however data type violations and specific values will be caught during the validation process and cause the model not to run (e.g set_diam: 0).


### Running the Model

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

#### Running from command line

1. Set current directory to **`BeRCM/model`**:

```bash
 cd /path/to/BeRCM/model
```

2. a. Run a _single_ instance of the model:

    ```bash
    python3 run.py PARAM_FILE
    ```
    Where `PARAM_FILE` is the path to the desired parameters file. To use the default parameters file set `PARAM_FILE=`**`parameters.yaml`**. 

    b. Run _multiple_ instances of the model:

    ```bash
    python3 parallel.py NUM_PROCESSES PARAM_FILE...
    ```
    Note that multiple parameter files can be passed to the **`parallel.py`** script. If more than one parameter file is passed then the number of files passed must be equal to the number of processes requested. 


#### Running in Spyder (THIS SECTION IS WIP)

1. Open **`run.py`** and **`parameters.py`** in Spyder

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

Ensuring `to_file` is set to `False` for one or both calls. Then follow steps 1 and 2 above.

### Results

All relevant information from the run will be stored in a Python [shelf](https://docs.python.org/3/library/shelve.html) in the **`plots`** directory:

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
    
    run-info-203012-2412-0142-1a4c536f-4bf3-4534-becf-90f207caf9e

For each run, a shelf file will contain the following information:

- Parameters

- Bed particle array

- For each iteration x:

    - Model particle array (at end of iteration x)
    
    - Available vertices (at end of iteration x)

    - Particles chosen for entrainment

- List of overall particle flux

- List of average particle age and range

### Plotting (THIS SECTION IS WIP)

**_Note_**: If doing another multiple runs, it is strongly suggested that you move the plots, rename them, or delete them. The model will overwrite the previous runs' plots.