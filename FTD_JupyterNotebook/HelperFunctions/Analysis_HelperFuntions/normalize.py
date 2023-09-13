import numpy as np

def normalize(arr, t_min, t_max):
    norm_arr = []

    diff = t_max - t_min
    diff_arr = np.nanmax(arr) - np.nanmin(arr)   

    for i in arr:
        temp = (((i - np.nanmin(arr))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
        
    return norm_arr