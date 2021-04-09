import sys
import os
import numpy as np
from scipy.stats import pareto
from scipy.spatial.distance import pdist, squareform



'''
Random stuff which is useful
'''



#Top level directory of code
headfolder  = os.path.dirname(os.path.realpath(__file__))
print(headfolder)
persloopfolder = None
assert persloopfolder is not None, "Please specify persloop location"


#Format and dpi of figures
fmt = 'png'
dpi = 300
fsize=None




#Make a path if it doesn't exist
def make_path(path):
    
    import errno
    
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


#Work out sensible default options of max_edge and min_pers if one/both is not specified
def default_maxedge(data):

    max_dist = pdist(data.T).max()
    max_edge = np.around(max_dist/2.,2)

    return max_edge


def default_minpers(data, name):

    scaling = 1. #default option

    #Dataset dependent scaling. Somewhat arbitrary: we want to make sure to only look for loops well above the gridscale
    if name == 'Lorenz':
        scaling = 10.

    if 'CDV' in name:
    #if name == 'CDV':
        #scaling = 130.
        scaling = 30.

    
    #Default is to first compute the pairwise distances, take the minimum
    #across all the pairwise distances (reducing the NxN matrix to N) and then
    #the mean across those (collapsing to a scalar)
    dist_arr = pdist(data.T)
    dist_arr_mat = squareform(dist_arr) + np.identity(data.shape[1])*dist_arr.max()
    min_dist = dist_arr_mat.min(axis=0).mean()
    min_pers = np.around(min_dist * scaling, 2)

    return min_pers




#Go from month string to number
def month2num(month):


    month2num_dict = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr' : 4, 'May' : 5, 'Jun' : 6,\
                      'Jul' : 7, 'Aug' : 8, 'Sep' : 9, 'Oct' : 10, 'Nov' : 11, 'Dec' : 12}

    return month2num_dict[month]

#Go from number to month string
def num2month(num):

    num2month_dict = {1 : 'Jan', 2 : 'Feb', 3 : 'Mar', 4 : 'Apr', 5 : 'May', 6 : 'Jun',\
                      7 : 'Jul', 8 : 'Aug', 9 : 'Sep', 10 : 'Oct', 11 : 'Nov', 12 : 'Dec'}

    return num2month_dict[num]


