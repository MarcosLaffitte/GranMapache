############################################################
#                                                          #
# - GranMapache: GRAphs-and-Networks MAPping Applications  #
#   with Cython and HEuristics                             #
#                                                          #
# - Description: suite for the analysis of bijections,     #
#   morphisms, alignments and other maps between graphs.   #
#                                                          #
# - Example of usage                                       #
#                                                          #
############################################################


# dependencies #############################################


# already in python ----------------------------------------
import time


# custom dependencies --------------------------------------
from .common_subgraphs import bla
from .integerization import ble
from .verbosity import bli


# main #####################################################


if __name__ == "__main__":
    # interpreted python
    initial_time = time.time()
    A = primes_python(1000)
    final_time = time.time()
    run_time = final_time - initial_time
    print("- Time (Python Interpreted): ", run_time, " seconds")
    # compiled python
    initial_time = time.time()
    A = primes_python_comp(1000)
    final_time = time.time()
    run_time = final_time - initial_time
    print("- Time (Python Compiled): ", run_time, " seconds")
    # cython
    initial_time = time.time()
    A = primes(1000)
    final_time = time.time()
    run_time = final_time - initial_time
    print("- Time (Cython): ", run_time, " seconds")
    # cython - submodule
    initial_time = time.time()
    A = primes_sub(1000)
    final_time = time.time()
    run_time = final_time - initial_time
    print("- Time (Cython Submodule): ", run_time, " seconds")


############################################################
############################################################
