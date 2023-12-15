from pycirclize import Circos
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fixedDensity as fd

def get_key_from_value(my_dict, search_value):
    for key, value in my_dict.items():
        if value == search_value:
            return key
    # If the value is not found, you can handle it as needed (e.g., return None or raise an exception)
    return None

def circularPlot(covMat, data_label, labels_dict, order, colordict, figTitle, fig_folder, 
                 long_track, short_track, region_spe = None,
                 fd_val = 10, 
                 cmp = False, 
                 fixed_density = False, 
                 long = False, short = False, region=False):
    
    legend_text = []
    ######################## Doing Group Comparison ########################
    if cmp:
        # Invert so that smaller p-value have greater thickness
        covMat_inv = 1 - covMat
        
        # Add ones to diagnoal values (this is because NaN values are in the diagnoal part)
        covMat_inv[np.identity(covMat_inv.shape[0], dtype=bool)] = 1
        
        # Substitute NaN values in covMat_inv to 0
        covMat_inv[np.isnan(covMat_inv)] = 0

        # Get upper triangle of covMat
        covMat_inv = np.triu(covMat_inv)

        # DataFrame of the Upper Triangle CovMat
        matrix_df = pd.DataFrame(covMat_inv, index=data_label, columns=data_label)

        # Check if Fixed Density is implemented
        if fixed_density:
            # Get the fixed density matrix with smallest fd_val
            fd_covMat = fd.fixedDensity(covMat, fd_val)

            # fixed density DataFrame
            df = pd.DataFrame(fd_covMat, index=data_label, columns=data_label)
            
            # get row, col indices for non-zero values
            row_indices, col_indices = np.where(df > 0)
            indices_list = list(zip(row_indices, col_indices))

            # convert this into Label from - to list
            label_list = []
            for f, t in indices_list:
                label_list.append((labels_dict[f], labels_dict[t]))

            if long: # Long Track + FD
                common_elements = [element for element in label_list if element in long_track]

                
                def link_kws_handler(from_label: str, to_label: str):
                    if (from_label, to_label) in common_elements:
                        # Set alpha, zorder values higher than other links for highlighting
                        return dict(alpha=1, zorder=1.0)
                    else:
                        return dict(alpha=0.1, zorder=0)

                # Initialize from matrix
                circos = Circos.initialize_from_matrix(
                    matrix_df,
                    space=3, # Degree between Nodes
                    r_lim=(93, 100), # Node width (size)
                    cmap = colordict, # Color Map
                    order = order,
                    label_kws=dict(r=100, size=12, color="black"), # Node Name
                    link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                    link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
                )

                # Print Original P-values
                legend_text = []
                for (from_label, to_label) in common_elements:
                    legend_text.append(f"{from_label}-{to_label}  p-value: {covMat[get_key_from_value(labels_dict, from_label)][get_key_from_value(labels_dict, to_label)]:.5f}")
                
            elif short: # Short Track + FD
                common_elements = [element for element in label_list if element in short_track]
               
                
                def link_kws_handler(from_label: str, to_label: str):
                    if (from_label, to_label) in common_elements:
                        # Set alpha, zorder values higher than other links for highlighting
                        return dict(alpha=1, zorder=1.0)
                    else:
                        return dict(alpha=0.1, zorder=0)

                # Initialize from matrix
                circos = Circos.initialize_from_matrix(
                    matrix_df,
                    space=3, # Degree between Nodes
                    r_lim=(93, 100), # Node width (size)
                    cmap = colordict, # Color Map
                    order = order,
                    label_kws=dict(r=100, size=12, color="black"), # Node Name
                    link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                    link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
                )

                # Print Original P-values
                legend_text = []
                for (from_label, to_label) in common_elements:
                    legend_text.append(f"{from_label}-{to_label}  p-value: {covMat[get_key_from_value(labels_dict, from_label)][get_key_from_value(labels_dict, to_label)]:.5f}")
                
            else: # Just FD
                # Define link_kws handler function to customize each link property
                def link_kws_handler(from_label: str, to_label: str):
                    if (from_label, to_label) in label_list:
                        # Set alpha, zorder values higher than other links for highlighting
                        return dict(alpha=1, zorder=1.0)
                    else:
                        return dict(alpha=0.1, zorder=0)

                # Initialize from matrix
                circos = Circos.initialize_from_matrix(
                    matrix_df,
                    space=3, # Degree between Nodes
                    r_lim=(93, 100), # Node width (size)
                    cmap = colordict, # Color Map
                    order = order,
                    label_kws=dict(r=100, size=12, color="black"), # Node Name
                    link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                    link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
                )

                # Print Original P-values
                legend_text = []
                for (from_ind, to_ind) in indices_list:
                    legend_text.append(f"{labels_dict[from_ind]}-{labels_dict[to_ind]} p-value: {covMat[from_ind][to_ind]:.5f}")
            
        # Fixed Density Not Implemented
        else:
            # Initialize from matrix
            circos = Circos.initialize_from_matrix(
                matrix_df,
                space=3, # Degree between Nodes
                r_lim=(93, 100), # Node width (size)
                cmap = colordict, # Color Map
                order = order,
                label_kws=dict(r=100, size=12, color="black"), # Node Name
                link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
            )
         
    ######################## JUST HC, Tau, TDP ########################
    else:
        covMat_norm = covMat

        # Substitute NaN values in covMat to 0
        covMat_norm[np.isnan(covMat_norm)] = 0

        # Get upper triangle of covMat
        covMat_norm = np.triu(covMat_norm)

        # DataFrame of the Upper Triangle CovMat
        matrix_df = pd.DataFrame(covMat_norm, index=data_label, columns=data_label)
            
        if long: # Long Track
            def link_kws_handler(from_label: str, to_label: str):
                if (from_label, to_label) in long_track:
                    # Set alpha, zorder values higher than other links for highlighting
                    return dict(alpha=1, zorder=1.0)
                else:
                    return dict(alpha=0.1, zorder=0)
                
            # Initialize from matrix
            circos = Circos.initialize_from_matrix(
                matrix_df,
                space=3, # Degree between Nodes
                r_lim=(93, 100), # Node width (size)
                cmap = colordict, # Color Map
                order = order,
                label_kws=dict(r=100, size=12, color="black"), # Node Name
                link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
            )
            
            # Print Covariance Values
            legend_text = []
            for (from_label, to_label) in long_track:
                legend_text.append(f"{from_label}-{to_label} Covariance Value: {covMat_norm[get_key_from_value(labels_dict, from_label)][get_key_from_value(labels_dict, to_label)]:.5f}")
                
        elif short: # Short Track
            def link_kws_handler(from_label: str, to_label: str):
                if (from_label, to_label) in short_track:
                    # Set alpha, zorder values higher than other links for highlighting
                    return dict(alpha=1, zorder=1.0)
                else:
                    return dict(alpha=0.1, zorder=0)
                
            # Initialize from matrix
            circos = Circos.initialize_from_matrix(
                matrix_df,
                space=3, # Degree between Nodes
                r_lim=(93, 100), # Node width (size)
                cmap = colordict, # Color Map
                order = order,
                label_kws=dict(r=100, size=12, color="black"), # Node Name
                link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
            )
            
            # Print Covariance Values
            legend_text = []
            for (from_label, to_label) in short_track:
                legend_text.append(f"{from_label}-{to_label} Covariance Value: {covMat_norm[get_key_from_value(labels_dict, from_label)][get_key_from_value(labels_dict, to_label)]:.5f}")
                
        else: # All
            if region: # Region specification
                def link_kws_handler(from_label: str, to_label: str):
                    if (from_label, to_label) in region_spe:
                        # Set alpha, zorder values higher than other links for highlighting
                        return dict(alpha=1, zorder=1.0)
                    else:
                        return dict(alpha=0.1, zorder=0)
                
                # Initialize from matrix
                circos = Circos.initialize_from_matrix(
                    matrix_df,
                    space=3, # Degree between Nodes
                    r_lim=(93, 100), # Node width (size)
                    cmap = colordict, # Color Map
                    order = order,
                    label_kws=dict(r=100, size=12, color="black"), # Node Name
                    link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                    link_kws_handler=link_kws_handler, # Highlighting Fixed Density Ones
                )

                # Print Covariance Values
                legend_text = []
                for (from_label, to_label) in region_spe:
                    legend_text.append(f"{from_label}-{to_label} Covariance Value: {covMat_norm[get_key_from_value(labels_dict, from_label)][get_key_from_value(labels_dict, to_label)]:.5f}")
                
            else:
                # Initialize from matrix
                circos = Circos.initialize_from_matrix(
                    matrix_df,
                    space=3, # Degree between Nodes
                    r_lim=(93, 100), # Node width (size)
                    cmap = colordict, # Color Map
                    order = order,
                    label_kws=dict(r=100, size=12, color="black"), # Node Name
                    link_kws=dict(direction=2, ec="black", lw=0.7), # Bidirectional, edgecolor, edge thickness
                )
    
    # Display Figure
    fig = circos.plotfig()

    # Title
    fig.suptitle(figTitle, fontsize=16)
    
    value_list = []
    for text in legend_text:
        val = float(text.split(':')[-1])
        value_list.append(val)
    
    value_list = np.array(value_list)
    
    # Add mean and std
    legend_text.append(f"Mean: {np.mean(value_list)}")
    legend_text.append(f"Std: {np.std(value_list)}")

#     plt.legend(legend_text, loc='upper right', handlelength=0, bbox_to_anchor=(1.5, 1.5))
    for i, line in enumerate(legend_text):
        plt.figtext(0.5, 0.01 - i * 0.02, line, ha='center', va='center')
    
    # Save Figure
    plt.savefig(fig_folder + '/' + figTitle + '.png')

    plt.show()