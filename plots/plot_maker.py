import numpy as np
import json 

import plotting


# Path to run-info file
fp_in = 'run-info-202105-3120-5743-3c8434c2-8ee4-4898-acae-9d1cb7f0e21a.json'
# Path to output directory
fp_out = './'
# Range to plot
iteration_range = [100, 150] 
# Plot height
y_limit = 10

with open(fp_in, 'r') as fp:
    run_info = json.load(fp)

for iter in range(iteration_range[0], iteration_range[1]+1):
    plotting.stream(iter, np.array(run_info["bed"]), np.array(run_info[str(iter)][0]), 
                                run_info["param"]["x_max"], y_limit, np.array(run_info[str(iter)][1]), fp_out)


