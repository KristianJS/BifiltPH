################################################################################
#                                                                              #
#        TOPOLOGICAL ANALYSIS OF DYNAMICAL SYSTEMS                             #
#                                                                              #
################################################################################


Code to carry out analysis of topological structure in the phase space of
a dynamical system. In particular, to compute persistent 1-cycles (i.e. loops)
and visualise sensible representatives of these. Done using a bifiltration
across density-distance space.

This code was the basis for all the analysis carried out in the paper
"A topological perspective on regimes in dynamical systems" by 
K. Strommen, M. Chantry, J. Dorrington and N. Otter (2020).


    *** BASICS ***

The main script is bifiltration.py, which loads a raw dataset, filters it by
density, computes a bifiltration of homology, and then generates basic plots.
All output and figures automatically stored. There are optargs that are helpful
when you run multiple times. Type

 > bifiltration.py -h

for explanations of arguments.


    *** HOW TO USE ***

* Install the conda environment using the supplied yaml file and make sure
  to activate it before running.

* Install PersLoop.

* In misc.py, change the headfolder variable to the full path of wherever
  you are putting this repo.

* Run > python bifiltration.py DATASETNAME --optargs


When adding a new dataset, you need to update loader.py appropriately.
It is also strongly recommended to add customised parameter choices in
bifiltration.py, since the default parameters may give bad performance
for your dataset. Similarly, in plots.py, parameter choices should be set,
which handle what angle to view the dataset from in 3D etc.

Please forgive idiosyncracies of this code. An attempt was made to balance
flexibility, ease of use and readability. No doubt we failed in many regards.


==================================================================================
Authors: K. Strommen, M. Chantry, J. Dorrington, N. Otter
==================================================================================
