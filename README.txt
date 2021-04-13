################################################################################
#                                                                              #
#            TOPOLOGICAL ANALYSIS OF DYNAMICAL SYSTEMS                         #
#                                                                              #
################################################################################


Code to carry out analysis of topological structure in the phase space of
a dynamical system. In particular, to compute persistent 0-cycles (i.e. components)
and persistent 1-cycles (i.e. loops) and visualise sensible representatives 
of these. Done using a bifiltration across density-distance space.

This code was the basis for all the analysis carried out in the paper
"A topological perspective on regimes in dynamical systems" by 
K. Strommen, M. Chantry, J. Dorrington and N. Otter (2021).


      *** BASICS ***

The main script is bifiltration.py, which loads a raw dataset (with loader.py), 
filters it by density (with filter.py), computes a bifiltration of homology
(with compute.py), and then generates basic plots (with plots.py).

All output and figures are automatically stored. There are optargs that are helpful
when you run multiple times. Type

 > bifiltration.py -h

for explanations of arguments, or read the scripts themselves, which are annotated.



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
  > cd bin
  > cmake ../perloop-src
  > cmake --build .
  
  NB: Persloop requires the boost libraries to install.


* In misc.py, you can change the headfolder variable to the full path of wherever
  you are putting this repo (it defaults to where the misc script is). You also 
  must update the persloopfolder variable to wherever you put the PersLoop executable,

  e.g. persloopfolder = "~/Persloop-viewer/bin/src/"

  Note it needs to point to the one in /bin/src to work.


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



     *** ACKNOWLEDGEMENTS ***

Special thanks to Hannah Christensen for allowing us to upload the Lorenz 96
data she generated, as described in the following papers:

* Arnold H. M., Moroz I. M. and Palmer T. N. (2013), 
  Stochastic parametrizations and model uncertainty in the Lorenz ’96 system,
  Phil. Trans. R. Soc. A.3712011047920110479,
   https://doi.org/10.1098/rsta.2011.0479

* Christensen, H.M., Moroz, I.M. & Palmer, T.N. (2015), 
  Simulating weather regimes: impact of stochastic and perturbed parameter schemes
  in a simple atmospheric model. Clim Dyn 44, 2195–2214,
  https://doi.org/10.1007/s00382-014-2239-9


==================================================================================
Authors: K. Strommen, M. Chantry, J. Dorrington, N. Otter
==================================================================================
