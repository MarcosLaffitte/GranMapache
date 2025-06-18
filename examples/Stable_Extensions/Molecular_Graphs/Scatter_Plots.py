################################################################################
#                                                                              #
#  - Analysis of Random Graphs for Extension of Partial Atonm-to-Atom Maps     #
#                                                                              #
#  - Plotting                                                                  #
#                                                                              #
#  - for WABI 2025                                                             #
#                                                                              #
#  - Made by Marcos Laffitte - Github @MarcosLaffitte                          #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import math
import time
import random
import pickle
from copy import deepcopy


# additional dependencies ------------------------------------------------------
import numpy as np
import cython
import networkx as nx
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


# parameters ###################################################################


# input ------------------------------------------------------------------------
time_gm_extender = None
time_gm_isomorphism = None
time_nx_isomorphism = None


# other variables --------------------------------------------------------------
X = []
C = []
D = []
time_gm_extension = []
time_gm_isomorphism = []
time_nx_isomorphism = []


# plot attributes --------------------------------------------------------------
plt.rcParams.update({"font.weight": "light", "font.family": "serif"})


# analysis #####################################################################



# task message
print("\n")
print("> Making Scatter Plots")
print("\n")



# open the file with graphs of each order
file_name = "data.csv"
with open(file_name, "r") as handle:
    times_by_reaction = handle.readlines()



# get number of nodes
for each_line in times_by_reaction[1:]:

    # split line
    values = each_line.split(",")

    # get parameters of graphs
    n_ITS = int(values[6])
    m_ITS = int(values[7])
    density_ITS = 100 * (m_ITS / (n_ITS * (n_ITS - 1) / 2))

    n_RC = int(values[8])
    m_RC = int(values[9])


    # get density
    X.append(100*n_RC/n_ITS)
    C.append(density_ITS)
    D.append(100*m_RC/m_ITS)


    # data for gm_extension
    T = float(values[2])*1000
    T = math.log(T)
    time_gm_extension.append(T)


    # data for gm_isomorphism
    T = float(values[4])*1000
    T = math.log(T)
    time_gm_isomorphism.append(T)


    # data for nx_isomorphism
    T = float(values[3])*1000
    T = math.log(T)
    time_nx_isomorphism.append(T)



# make scatter plot gm_extension
plt.scatter(X, time_gm_extension, c = C, alpha = 0.5, cmap = cmr.get_sub_cmap("viridis", 0.2, 0.8))

# figure attributes
plt.xlabel("Proportion of reacting vertices [%]", labelpad = 10, size = 20)
plt.ylabel("Log of running time [ms]", labelpad = 10, size = 20)
plt.ylim(-0.8, 3.8)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
plt.tight_layout()

# save figure
file_name = "Scatter_gm_extension.png"
plt.savefig(file_name)
plt.close()



# make scatter plot nx_isomorphism
plt.scatter(X, time_nx_isomorphism, c = C, alpha = 0.5, cmap = cmr.get_sub_cmap("viridis", 0.2, 0.8))

# figure attributes
plt.xlabel("Proportion of reacting vertices [%]", color = "w", labelpad = 10, size = 20)
plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 20)
plt.ylim(-0.8, 3.8)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
plt.tight_layout()

# save figure
file_name = "Scatter_nx_isomorphism.png"
plt.savefig(file_name)
plt.close()



# make scatter plot gm_isomorphism
scp = plt.scatter(X, time_gm_isomorphism, c = C, alpha = 0.5, cmap = cmr.get_sub_cmap("viridis", 0.2, 0.8))

# figure attributes
plt.xlabel("Proportion of reacting vertices [%]", color = "w", labelpad = 10, size = 20)
plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 20)
plt.ylim(-0.8, 3.8)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)

# save figure
plt.tight_layout()
file_name = "Scatter_gm_isomorphism.png"
plt.savefig(file_name)
# plt.close()



# plot color bar of scatter plots
cbar = plt.colorbar(scp)
cbar.solids.set(alpha = 1)
cbar.ax.set_ylabel("Proportion of edges in ITS [%]", labelpad = 25, rotation = 270, size = 18)
cbar.ax.tick_params(labelsize = 16)
plt.tight_layout()

# save figure
file_name = "Scatter_gm_isomorphism_with_colorbar.png"
plt.savefig(file_name)
plt.close()



# line jump message
print("> Finished")
print("\n")


################################################################################
################################################################################
