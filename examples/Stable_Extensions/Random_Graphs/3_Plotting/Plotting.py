################################################################################
#                                                                              #
#  - Analysis of Random Graphs for Extension of Partial Atonm-to-Atom Maps     #
#                                                                              #
#  - Plotting                                                                  #
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
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


# parameters ###################################################################


# input ------------------------------------------------------------------------
in_time_gm_extender = None
in_time_gm_isomorphism = None
in_time_nx_isomorphism = None
in_time_nx_ismags_algo = None


# parameters of remainder graphs -----------------------------------------------
nodes_remainders = [100, 125, 150, 175, 200]
# NOTE: density of at least 2 is needed for connectedness with 100 nodes,
# but with density exactly 2 and 125 nodes it seems way too hard to find a random connected graph,
# and density of at most 97 is also needed for possibility of construction with 100 vertices,
# thus the 3-97 range for density was selected, which seems to work well
densities_remainders = list(range(3, 97 + 1))
density_counter = 0


# other variables --------------------------------------------------------------
window_size = 5
max_rolling = 0
window_gm_extender = []
window_gm_isomorphism = []
window_nx_isomorphism = []
window_nx_ismags_algo = []
lower_window_gm_extender = []
lower_window_gm_isomorphism = []
lower_window_nx_isomorphism = []
lower_window_nx_ismags_algo = []
upper_window_gm_extender = []
upper_window_gm_isomorphism = []
upper_window_nx_isomorphism = []
upper_window_nx_ismags_algo = []
times_by_order = []
slice_gm_extender = []
slice_gm_isomorphism = []
slice_nx_isomorphism = []
slice_nx_ismags_algo = []
each_list_of_times_for_fixed_density = []


# plot data --------------------------------------------------------------------
average_times_by_order_gm_extender = []
average_times_by_order_gm_isomorphism = []
average_times_by_order_nx_isomorphism = []
average_times_by_order_nx_ismags_algo = []
stddev_times_by_order_gm_extender = []
stddev_times_by_order_gm_isomorphism = []
stddev_times_by_order_nx_isomorphism = []
stddev_times_by_order_nx_ismags_algo = []
minimum_times_by_order_gm_extender = []
minimum_times_by_order_gm_isomorphism = []
minimum_times_by_order_nx_isomorphism = []
minimum_times_by_order_nx_ismags_algo = []
maximum_times_by_order_gm_extender = []
maximum_times_by_order_gm_isomorphism = []
maximum_times_by_order_nx_isomorphism = []
maximum_times_by_order_nx_ismags_algo = []
rolling_average_times_by_order_gm_extender = []
rolling_average_times_by_order_gm_isomorphism = []
rolling_average_times_by_order_nx_isomorphism = []
rolling_average_times_by_order_nx_ismags_algo = []


# plot attributes --------------------------------------------------------------
plt.rcParams.update({"font.weight": "light", "font.family": "serif"})
# ylims = {100: (1.4, 5.2),
#          125: (1.4, 5.7),
#          150: (2.0, 6.2),
#          175: (2.3, 6.5),
#          200: (2.5, 7.0)}
ylims = {100: (1.4, 9.5),
         125: (1.4, 9.5),
         150: (2.0, 9.5),
         175: (2.3, 9.5),
         200: (2.5, 9.5)}


# analysis #####################################################################


# task message
print("\n")
print("> Plotting times")
print("\n")


# get number of nodes
for each_order in nodes_remainders:

    # open the file with graphs of each order
    file_name = "Times_Random_Graphs_" + str(each_order) + "_nodes.pkl"
    with open(file_name, "rb") as handle:
        times_by_order = pickle.load(handle)

    # reinitialize variables
    average_times_by_order_gm_extender = []
    average_times_by_order_gm_isomorphism = []
    average_times_by_order_nx_isomorphism = []
    average_times_by_order_nx_ismags_algo = []
    stddev_times_by_order_gm_extender = []
    stddev_times_by_order_gm_isomorphism = []
    stddev_times_by_order_nx_isomorphism = []
    stddev_times_by_order_nx_ismags_algo = []

    # open each sublist of graphs by density
    for each_list_of_times_for_fixed_density in times_by_order:

        # average for gm_extender
        slice_gm_extender = [math.log(a*1000) for (a, b, c, d) in each_list_of_times_for_fixed_density]
        average_times_by_order_gm_extender.append(sum(slice_gm_extender)/len(slice_gm_extender))
        stddev_times_by_order_gm_extender.append(np.std(slice_gm_extender))
        minimum_times_by_order_gm_extender.append(min(slice_gm_extender))
        maximum_times_by_order_gm_extender.append(max(slice_gm_extender))

        # average for gm_isomorphism
        slice_gm_isomorphism = [math.log(b*1000) for (a, b, c, d) in each_list_of_times_for_fixed_density]
        average_times_by_order_gm_isomorphism.append(sum(slice_gm_isomorphism)/len(slice_gm_isomorphism))
        stddev_times_by_order_gm_isomorphism.append(np.std(slice_gm_isomorphism))
        minimum_times_by_order_gm_isomorphism.append(min(slice_gm_isomorphism))
        maximum_times_by_order_gm_isomorphism.append(max(slice_gm_isomorphism))

        # average for nx_isomorphism
        slice_nx_isomorphism = [math.log(c*1000) for (a, b, c, d) in each_list_of_times_for_fixed_density]
        average_times_by_order_nx_isomorphism.append(sum(slice_nx_isomorphism)/len(slice_nx_isomorphism))
        stddev_times_by_order_nx_isomorphism.append(np.std(slice_nx_isomorphism))
        minimum_times_by_order_nx_isomorphism.append(min(slice_nx_isomorphism))
        maximum_times_by_order_nx_isomorphism.append(max(slice_nx_isomorphism))

        # average for nx_ismags_algo
        slice_nx_ismags_algo = [math.log(d*1000) for (a, b, c, d) in each_list_of_times_for_fixed_density]
        average_times_by_order_nx_ismags_algo.append(sum(slice_nx_ismags_algo)/len(slice_nx_ismags_algo))
        stddev_times_by_order_nx_ismags_algo.append(np.std(slice_nx_ismags_algo))
        minimum_times_by_order_nx_ismags_algo.append(min(slice_nx_ismags_algo))
        maximum_times_by_order_nx_ismags_algo.append(max(slice_nx_ismags_algo))

    # plot in density
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.plot(densities_remainders, average_times_by_order_gm_extender, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_extender")
    plt.plot(densities_remainders, average_times_by_order_gm_isomorphism, linewidth = 1,  marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.plot(densities_remainders, average_times_by_order_nx_isomorphism, linewidth = 1,  marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.plot(densities_remainders, average_times_by_order_nx_ismags_algo, linewidth = 1,  marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left", prop = font, fontsize = 15)

    # figure attributes
    plt.title(str(each_order) + " nodes in remainder graph", fontsize = 12, fontweight = "bold", fontfamily = "serif")
    if(each_order in [100, 125]):
        plt.xlabel("Proportion of edges in remainder graph [%]", labelpad = 10, size = 18)
        plt.ylabel("Log of running time [ms]", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Average_Times_Random_Graphs_" + str(each_order) + "_nodes.pdf"
    plt.savefig(file_name)
    plt.close()

    # reinitialize variables
    rolling_average_times_by_order_gm_extender = []
    rolling_average_times_by_order_gm_isomorphism = []
    rolling_average_times_by_order_nx_isomorphism = []
    rolling_average_times_by_order_nx_ismags_algo = []

    # get rolling average
    max_rolling = len(times_by_order)
    for k in range(max_rolling):

        # get lower window
        if(k == 0):
            lower_window_gm_extender = []
            lower_window_gm_isomorphism = []
            lower_window_nx_isomorphism = []
            lower_window_nx_ismags_algo = []
        else:
            lower_window_gm_extender = average_times_by_order_gm_extender[max(0, k - window_size) : k]
            lower_window_gm_isomorphism = average_times_by_order_gm_isomorphism[max(0, k - window_size) : k]
            lower_window_nx_isomorphism = average_times_by_order_nx_isomorphism[max(0, k - window_size) : k]
            lower_window_nx_ismags_algo = average_times_by_order_nx_ismags_algo[max(0, k - window_size) : k]

        # get upper window
        if(k == max_rolling - 1):
            upper_window_gm_extender = []
            upper_window_gm_isomorphism = []
            upper_window_nx_isomorphism = []
            upper_window_nx_ismags_algo = []
        else:
            upper_window_gm_extender = average_times_by_order_gm_extender[k+1 : min(k + window_size, max_rolling - 1) + 1]
            upper_window_gm_isomorphism = average_times_by_order_gm_isomorphism[k+1 : min(k + window_size, max_rolling - 1) + 1]
            upper_window_nx_isomorphism = average_times_by_order_nx_isomorphism[k+1 : min(k + window_size, max_rolling - 1) + 1]
            upper_window_nx_ismags_algo = average_times_by_order_nx_ismags_algo[k+1 : min(k + window_size, max_rolling - 1) + 1]

        # get full windows
        window_gm_extender = lower_window_gm_extender + [average_times_by_order_gm_extender[k]] + upper_window_gm_extender
        window_gm_isomorphism = lower_window_gm_isomorphism + [average_times_by_order_gm_isomorphism[k]] + upper_window_gm_isomorphism
        window_nx_isomorphism = lower_window_nx_isomorphism + [average_times_by_order_nx_isomorphism[k]] + upper_window_nx_isomorphism
        window_nx_ismags_algo = lower_window_nx_ismags_algo + [average_times_by_order_nx_ismags_algo[k]] + upper_window_nx_ismags_algo

        # get averages
        rolling_average_times_by_order_gm_extender.append(sum(window_gm_extender) / len(window_gm_extender))
        rolling_average_times_by_order_gm_isomorphism.append(sum(window_gm_isomorphism) / len(window_gm_isomorphism))
        rolling_average_times_by_order_nx_isomorphism.append(sum(window_nx_isomorphism) / len(window_nx_isomorphism))
        rolling_average_times_by_order_nx_ismags_algo.append(sum(window_nx_ismags_algo) / len(window_nx_ismags_algo))






    # plot rolling average gm_extender
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_extender, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_extender")
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_ismags_algo, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left")

    # figure attributes
    plt.xlabel("Proportion of edges in ITS [%]", color = "w", labelpad = 10, size = 18)
    plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Rolling_Average_Times_Random_Graphs_" + str(each_order) + "_nodes_global.pdf"
    plt.savefig(file_name)
    plt.close()





    # plot rolling average gm_extender
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.errorbar(densities_remainders, rolling_average_times_by_order_gm_extender, linewidth = 1,
                 yerr = stddev_times_by_order_gm_extender, elinewidth = 0.5,
                 marker = ".", markersize = 2.5, label = "gm_extender")
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_ismags_algo, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left")

    # figure attributes
    plt.xlabel("Proportion of edges in ITS [%]", color = "k", labelpad = 10, size = 18)
    plt.ylabel("Log of running time [ms]", color = "k", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Rolling_Average_Times_Random_Graphs_" + str(each_order) + "_nodes_gm_extender.pdf"
    plt.savefig(file_name)
    plt.close()






    # plot rolling average gm_isomorphism
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_extender, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_extender")
    plt.errorbar(densities_remainders, rolling_average_times_by_order_gm_isomorphism, linewidth = 1,
                 yerr = stddev_times_by_order_gm_isomorphism, elinewidth = 0.5,
                 marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_ismags_algo, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left")

    # figure attributes
    plt.xlabel("Proportion of edges in ITS [%]", color = "w", labelpad = 10, size = 18)
    plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Rolling_Average_Times_Random_Graphs_" + str(each_order) + "_nodes_gm_isomorphism.pdf"
    plt.savefig(file_name)
    plt.close()






    # plot rolling average nx_isomorphism
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_extender, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_extender")
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.errorbar(densities_remainders, rolling_average_times_by_order_nx_isomorphism, linewidth = 1,
                 yerr = stddev_times_by_order_nx_isomorphism, elinewidth = 0.5,
                 marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_ismags_algo, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left")

    # figure attributes
    plt.xlabel("Proportion of edges in ITS [%]", color = "w", labelpad = 10, size = 18)
    plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Rolling_Average_Times_Random_Graphs_" + str(each_order) + "_nodes_nx_isomorphism.pdf"
    plt.savefig(file_name)
    plt.close()






    # plot rolling average nx_ismags_algo
    plt.text(3 - 1, ylims[each_order][1] + 0.05, "3%", fontsize = 12)
    plt.text(97 - 2, ylims[each_order][1] + 0.05, "97%", fontsize = 12)
    plt.vlines(3, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.vlines(97, ymin = ylims[each_order][0], ymax = ylims[each_order][1], color = "k", linestyle = "--", linewidth = 1.2)
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_extender, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_extender")
    plt.plot(densities_remainders, rolling_average_times_by_order_gm_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "gm_isomorphism")
    plt.plot(densities_remainders, rolling_average_times_by_order_nx_isomorphism, linewidth = 1, marker = ".", markersize = 2.5, label = "nx_isomorphism")
    plt.errorbar(densities_remainders, rolling_average_times_by_order_nx_ismags_algo, linewidth = 1,
                 yerr = stddev_times_by_order_nx_ismags_algo, elinewidth = 0.5,
                 marker = ".", markersize = 2.5, label = "nx_ismags_algo")

    # legend
    # plt.legend(loc = "upper left")

    # figure attributes
    plt.xlabel("Proportion of edges in ITS [%]", color = "w", labelpad = 10, size = 18)
    plt.ylabel("Log of running time [ms]", color = "w", labelpad = 10, size = 18)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(ylims[each_order][0], ylims[each_order][1])
    plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
    plt.tight_layout()

    # save figure
    file_name = "Rolling_Average_Times_Random_Graphs_" + str(each_order) + "_nodes_nx_ismags_algo.pdf"
    plt.savefig(file_name)
    plt.close()






# line jump message
print("> Finished")
print("\n")


################################################################################
################################################################################
