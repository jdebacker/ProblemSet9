#update r path after updating BQ path
# use market clearing
L_path = np.ones(T) * agg.get_L(n, omega_SS)
# find implied r
new_rpath = firm.get_r(L_path, new_Kpath, r_params)
# update r
r_path_init[:T] = xi * new_rpath[:T] + (1 - xi) * rpath_init[:T]

# change dist to the deviations for the BQ path and the r path.
dist = np.absolute(BQpath_init[:T] - new_BQpath[:T]).max()+np.absolute(rpath_init[:T] - new_rpath[:T]).max()

