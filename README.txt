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

* Install the conda environment using the supplied yaml file:
  
  > conda env create -f env.yaml

  Its default name is 'pholo'. Make sure to activate it before running:

  > conda activate pholo

* Install PersLoop using the following commands:

  > git clone https://github.com/Sayan-m90/Persloop-viewer.git
  > cd Persloop-viewer
  
  #Checkout the commit used in this work
  > git checkout 165b1412d0ec70128cc393477e2e1480708bbeed
  > mkdir bin
  > cmake ../perloop-src
  > cmake --build .

* In misc.py, change the headfolder variable to the full path of wherever
  you are putting this repo, and also update the persloopfolder variable
  to wherever you put PersLoop,

  e.g. persloopfolder = "~/Persloop-viewer/bin/src/"

* To run a generic example using existing data, you can do

  > python bifiltration.py Lorenz63

  For other optional arguments, see the help statement:

  > ./bifiltration.py -h



    *** ADDING NEW DATASETS ***


When adding a new dataset, you need to update loader.py appropriately.
It is also strongly recommended to add customised parameter choices in
bifiltration.py, since the default parameters may give bad performance
for your dataset. Similarly, in plots.py, parameter choices should be set,
which handle what angle to view the dataset from in 3D etc.



    *** CONTACT DETAILS AND BUGS ***

Please report bugs as an issue, or email kristianjstr@gmail.com

Please forgive idiosyncracies of this code. An attempt was made to balance
flexibility, ease of use and readability. No doubt we failed in many regards.


==================================================================================
Authors: K. Strommen, M. Chantry, J. Dorrington, N. Otter
==================================================================================
