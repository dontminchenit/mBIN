import numpy as np

def normalize2d(arr, t_min, t_max):
    norm_arr = np.zeros(arr.shape)

    # Difference between our set t_max and t_min
    diff = t_max - t_min
    # Difference between max and min value of our 2d ndarray
    arr_min = np.nanmin(arr)
    arr_max = np.nanmax(arr)
    diff_arr = arr_max - arr_min   

    for i in range(arr.shape[0]): # Row
        for j in range(arr.shape[1]): # Column
            temp = (((arr[i,j] - arr_min) * diff) / diff_arr) + t_min
            norm_arr[i, j] = temp
        
    return norm_arr
