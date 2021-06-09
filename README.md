# Py_BeRCM

### Requirements

Installing Python requirements:

```bash
 pip3 install -r requirements.txt
```

### Important Terminology in Model

- **Bed Particle**: A particle embedded in the stream bed. Does not move.
- **Model Particle**: A particle subject to entrainment events. Moves.
- **Event Particle**: A Model Particle which has been selected for entrainment.
- **Desired Hop**: The hop distance an Event Particle is assigned at the beginning of an entrainment event
- **Stream**: Consists of the Bed and Model particles.

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

2. Run the model:
    - Run the model _with_ stream plots:

    ```bash
    python3 run.py 
    ```

   - Run the model _without_ stream plots (Flux plots will still be created):

    ```bash
    python3 run.py --no-plots
    ```

#### Running in Spyder

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

### Setting Parameters

The model comes with a parameters file, `parameters.py`, which users can edit. The file is located in the **`model`** directory:

```
project
│   README.md
│   requirements.txt
└───model
│   │   parameters.py
│   │   run.py
│   │   ...
```

Comments indicate what values have been tested for each paramter and within what range. Users are welcome to enter parameters beyond and below these bounds to see the results, however some values will be caught and the model will not run (i.e grain diameter < 0).

### Results

You can view the plots in the **`plots`** directory:

```
RCM_BedloadTransport
│   README.md
│   requirements.txt
└───model
│   │   ...
└───plots <--- here!
│   │   ...
```

If you chose to run the model with plotting, then this directory will contain plots of the stream during each iteration as well as the final flux histogram. If you ran the model without plotting then the flux histogram is the only plot produced.

**_Note_**: If doing another multiple runs, it is strongly suggested that you move the plots, rename them, or delete them. The model will overwrite the previous runs' plots.