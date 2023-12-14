import sys, os
import pandas as pd
import numpy  as np
from scipy.optimize import linear_sum_assignment

sim_file, num_file = os.path.join(sys.argv[1], 'sim.xls'), os.path.join(sys.argv[1], 'num.xls')
cost, cell_number_to_spot = pd.read_csv(sim_file).T.values, pd.read_csv(num_file)
cell_number_to_spot = [xx[0] for xx in cell_number_to_spot.values.tolist()]

location_repeat = np.zeros(cost.shape[1])
location_repeat = np.repeat(np.arange(len(cell_number_to_spot)), cell_number_to_spot)

# Assigning individual cells to spatial coordinates.
location_repeat = location_repeat.astype(int)
cost_mat = cost[location_repeat, :]
row_ind, y = linear_sum_assignment(cost_mat.T)	

assigned_nodes = location_repeat[y]
np.savetxt(os.path.join(sys.argv[1], 'output.txt'), assigned_nodes, fmt = '%d', delimiter = ',')
