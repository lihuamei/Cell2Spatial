import sys, os
import pandas as pd
import numpy  as np
from lapjv import lapjv

if __name__ == '__main__':
    tmp_dir, sim_file, num_file, prefix  = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    cost, cell_number_to_spot = pd.read_csv(sim_file).T.values, pd.read_csv(num_file)
    cell_number_to_spot = [xx[0] for xx in cell_number_to_spot.values.tolist()]

    location_repeat = np.zeros(cost.shape[1])
    location_repeat = np.repeat(np.arange(len(cell_number_to_spot)), cell_number_to_spot)

    # Assigning individual cells to spatial coordinates.
    location_repeat = location_repeat.astype(int)
    cost_mat = cost[location_repeat, :]

    np.random.seed(1234)
    cost_scaled = cost_mat +  1e-16 * np.random.rand(cost_mat.shape[0], cost_mat.shape[1])
    _, y, _ = lapjv(cost_scaled)

    assigned_nodes = location_repeat[y]
    np.savetxt(os.path.join(sys.argv[1], prefix + '_output.txt'), assigned_nodes, fmt = '%d', delimiter = ',')

