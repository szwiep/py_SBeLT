import numpy as np
import json 

import plotting


# Path to run-info file
fp_in = 'run-info-202106-0417-1051-db2284e0-493b-471a-b4c3-204a40598853.json'
# Path to output directory
fp_out = './hd_test2/'
# Range to plot
iteration_range = [99000, 99010] 
# Plot height
y_limit = 10

with open(fp_in, 'r') as fp:
    run_info = json.load(fp)

for iter in range(iteration_range[0], iteration_range[1]+1):
    plotting.stream(iter, np.array(run_info["bed"]), np.array(run_info[str(iter)][0]), 
                                run_info["param"]["x_max"], y_limit, np.array(run_info[str(iter)][1]), fp_out)


plotting.flux_info(np.array(run_info["flux"]), run_info["param"]["n_iterations"], fp_out)