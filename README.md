# Gran Mapache - welcome to the repo

<p align="center">
<strong>GranMapache:</strong> <strong>GRA</strong>phs-and-<strong>N</strong>etworks <strong>MAP</strong>ping <strong>A</strong>pplications with <strong>C</strong>ython and <strong>HE</strong>uristics
</p>

<p align="center">
Suite for the Analysis of Bijections, Morphisms,<br/>
Alignments and other Maps between Graphs
</p>


## ToDo / Working

> done: integerization

> done: good-bad map with maximum extension

> done: add wheel for packages

> done: statically type python objects

> done: statically type counters and bounds in for loops

> done: migrating partial_maps into C++ objects

> done: simplifying code with C++ structs

> done: implemented C++ unordered maps for fast (constant look-ups).

> done: implementing and testing iterative version of the algorithm(s) against current recursive versions.

> done: transition into unordered sets and unordered maps for scalability with constat look ups

> done: isomorphism

- working: simplify feasability testing

- working: add examples with isomorphisms to standalone

- working: mcs routines

- working: test isomorphism and mcs routines

- modules: integerization, partial_maps, common_subgraphs, alignment, isomorphism, rebalancing, verbosity, visualization, (maybe) classes, (maybe) tests.

- versions: pip for linux, pip for windows, stand-alone.


## Institutions

> Center for Scalable Data Analytics and Artificial Intelligence, Leipzig / Dresden, Germany. See <a href="https://scads.ai/">ScaDS.AI</a>.<br/>

> Bioinformatics Group, Department of Computer Science, Leipzig University, Germany. See <a href="https://www.bioinf.uni-leipzig.de/">Bioinf</a>.<br/>

## Developed by

- Marcos E. González Laffitte<br/>
  Bioinformatics Group and ScaDS.AI, Leipzig University<br/>
  marcoslaffitte@gmail.com<br/>
  marcos@bioinf.uni-leipzig.de<br/>

- Prof. Dr. Peter F. Stadler<br/>
  Bioinformatics Group and ScaDS.AI and Interdisciplinary Center for Bioinformatics, Leipzig University<br/>
  studla@bioinf.uni-leipzig.de<br/>


## Description

<div align="justify">
the description
</div>
<br/>

**Note:** a note.

## Cite as

This repository was developed as part of the contribution:

**[1]** TBA
> **Link:** https:

<div align="justify">
There you can find detailed information on the algorithms implemented here. This work was developed for research purposes. Please cite as above if you find this work or these programs useful for your own research.
</div>
<br/>

<div align="center">
<strong>Some figure from <a href="link">[1]</a></strong><br/>
</div>



## Instructions

###### 1) installing required python version into anaconda environment
```
conda create -n [env_name] python=3.11.10
```
###### 2) for all the following, activate anaconda environment with
```
conda activate [env_name]
```
###### 3) In order to install this package you need to have C and C++ compilers. If you do not have these already in your system you can install them INSIDE the anaconda environment with the following command:
```
conda install -c conda-forge cxx-compiler
```
###### 4) installing the package and dependencies INSIDE anaconda environment, go to the folder where setup.py is located and run
```
pip install .
```
###### 5) there is a python script with an example on how to call the functions in the packages, called Standalone_Example.py inside the standalone directory. To run this use:
```
python Standalone_Example.py
```
###### For uninstalling gmapache from anaconda environment
```
pip uninstall gmapache
```
###### For removing anaconda environmnet completely
```
conda remove -n [env_name] --all
```



## Additional Information

Some additional information


## LICENSE

The programs in this repository are part of the work in<br/>
[work]<br/>
and are released under<br/>
<strong>MIT License Copyright (c) 2024 Marcos E. González Laffitte</strong><br/>
See <a href="./LICENSE">LICENSE</a> file for full license details.