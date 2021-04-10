
#NOTES
#To do more than 3D pers you need to change threeD=False
#and compile the all-dim version of persloop
#this is unsupport by this work.

import gudhi as gd
import numpy as np
import sys
from misc import persloopfolder

peps = 10**(-7)
threeD = True
assert threeD  #PersLoop is not reliable or stable for dimensions>3


#Compute Euclidean distance between two points defined by data indices
def dist(ind1, ind2, data):

    d = np.linalg.norm(data[:,ind1]-data[:,ind2])
    return d


#Sort a number of generators returned from Gudhi into distinct components
def get_components(generators):
    l = generators
    out = []
    while len(l)>0:
        first, *rest = l
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        out.append(first)
        l = rest

    points = list(out)

    h0_arr = []
    for n in range(len(points)):
        points_l = list(points[n])
        points_l.sort()
        h0_arr.append(np.array(points_l))

    return h0_arr


#Get the N longest lived components of some dataset
#Input is a simplex tree object (filt) and its corresponding barcode
def oldest_components(data, filt, barcode, N):
    
    print(100*'-')
    print('       CONNECTED COMPONENTS:')
    
    #Get number of time indices
    datalength = data.shape[-1]
    print(" Datalength = %s" % datalength)

    filt_vals = np.array([filt[n][1] for n in range(len(filt))])
    death_times = []
    for component in barcode:
        if component[0] == 0: #we only want deathtimes of components (=H0, hence the 0 value)
            deathtime = component[1][1]
            death_times.append(deathtime)
       

    #If there are <N components to pick from then just take all of them
    if (len(death_times)<N):
        print(" Less components available than what was asked for. Redefining N.")
        N = len(death_times)
    
    print(" N = %s" % N)
    
    #Now add 0 to the end for convenience
    death_times.append(0.)
      
    #Get the generators up to this filtration value
    filt1 = death_times[N-1]
    #filt0 = death_times[N]
    filt0 = 0.
    print(' Looking at filtration values between %s and %s' % (filt1, filt0))
    generators_between = []
    for n in range(len(filt)):
        if (filt1>filt_vals[n]>filt0) and (len(filt[n][0]) > 1):
            generators_between.append(filt[n][0])
    
   

    #Get the first guess of connected components
    print(" ...building initial set of components from Gudhi generators...")
    components = get_components(generators_between)
   
    
    #Now add in all points within a ball of radius death_time around each generator of these components,
    #and each point within a ball of radius death_time around each of these new points, etc.,
    #recursively until you've exhausted all points in your dataset
    maximum_radius = death_times[N]
    
    missing_points = []
    for n in range(data.shape[-1]):
        found = False
        for comp in components:
            if n in comp:
                found = True
                
        if not(found):
            missing_points.append(n)
    
    missing_points.sort()
    M = len(missing_points)

    print(" Missing points: %s" % M)
    
    if M == 0:
        print(" All components found!")
    
        #We also want the number of points of each component
        comp_sizes = []
        for comp in components:
            comp_sizes.append(len(comp))
        comp_sizes = np.array(comp_sizes)
        print(" Size of components found: %s" % (comp_sizes))
        return components, death_times, comp_sizes
    else:
        pass
    
    print(" ...searching for singletons...")
    
    #First find singletons. These are points that are >maximum_radius radius away from any other point
    inds_to_pop = []
    for m in missing_points:
        dists = [dist(m, k, data) for k in range(data.shape[-1])]
        dists.pop(dists.index(0.))
        if np.array(dists).min() >= maximum_radius:
            print(" Found singleton with nearest neighbour at distance %s" % np.array(dists).min())
            components.append(np.array([m]))
            inds_to_pop.append(missing_points.index(m))
    
    #missing_points.pop(inds_to_pop)
    if len(inds_to_pop) > 0:
        for index in sorted(inds_to_pop, reverse=True):
            del missing_points[index]
    
    
    M = len(missing_points)
    
    def find_comp(m, components):
        
        for n in range(len(components)):
            comp = components[n]
            if m in comp:
                pass
            else:
                for ind in comp:
                    if dist(m, ind, data) <= maximum_radius:
                        return n


    cnter=0 #Don't loop more than 20 times
    while (M>0 and cnter < 20):
        print(" Still have missing values. Probing balls of increasing radii")
        print(" ...missing points: %s" % M)
        for m in missing_points:
            n = find_comp(m, components)
            if not(n is None) and not(m in components[n]):
                components[n] = np.append(components[n], m)
                missing_points.pop(missing_points.index(m))
                M -= 1
            cnter += 1
   
    if cnter == 20:
        print(" WARNING: Counter reached its limit!!")
        print(" Proceeding anyway...")
    else:
        print(" All components found!")


    #Kill any doubles
    components_clean = []
    for comp in components:
        components_clean.append(np.array(list(set(comp))))
    
    components_clean = components
    
    #We also want the number of points of each component
    comp_sizes = []
    for comp in components:
        comp_sizes.append(len(comp))
    comp_sizes = np.array(comp_sizes)

    print(" Size of components found: %s" % (comp_sizes))

    return components_clean, death_times, comp_sizes



def persloop(data_location=None,data=None,use_pl=True,**kwargs):
    #Desired arguments:
    # data or data_location
    # max_edge: based on computational resources
    # min_pers: distance from diag
    # sparse: recommend 0.8-0.9 for large datasets
    assert (data_location is not None or data is not None),'No data provided!'
    if 'nme' not in kwargs.keys():
        print("No name specified, using 'tmp'")
        kwargs['nme']='tmp'
    if data is None:
        #Do some data loading
        print("Need to write code for data loading! QUITTING")
        sys.exit()
        data = data
    proc4persloop(data,**kwargs)
    if use_pl:
        runpersloop(kwargs['nme'])
    return
        

def proc4persloop(data,nme,num_comp_to_keep=10,max_edge=10,min_pers=1,sparse=None,pre_sparse=None,tkbackend=False):
    assert int(gd.__version__[0]) >= 3
    #Uses Gudhi to process data and save data plus filtration
    if pre_sparse is not None:
        data = np.array(gd.sparsify_point_set(data,min_squared_dist=pre_sparse**2))
    print(100*'-')
    print(' *** STARTING HOMOLOGICAL COMPUTATIONS ***')

    #Immediately check datalength and cut it down to a maximum size of 20000
    maxlength = 20000
    datalength = data.shape[0]
    print(" Datalength = %s" % datalength)
    if (datalength > maxlength):
        print(" Data is too big! Subselecting first %s points" % maxlength)
        data = data[:maxlength,:]
    
    rc = gd.RipsComplex(data,max_edge_length=max_edge,sparse=sparse)
    Rips_simplex_tree = rc.create_simplex_tree(max_dimension=2)
    print(' Tree created')

    BarCodes_Rips0 = Rips_simplex_tree.persistence(min_persistence=min_pers)
    #print(BarCodes_Rips0)
    print(' Barcodes done')
    print('   Printing loops:')
    loop_births, loop_deaths = printloops(BarCodes_Rips0)
    Filt = Rips_simplex_tree.get_filtration()
    print(' Filtration gotten')

   
    print(' Making generic birth/death plot')
    import matplotlib
    if tkbackend:
        matplotlib.use('QT4Agg') #TKAgg
    import matplotlib.pyplot as plt
    gd.plot_persistence_diagram(BarCodes_Rips0)
    #plt.show()
    fname = nme + 'birthdeath.png' 
    plt.savefig(fname,dpi=400)
    print(" A birth/death plot has been saved to: %s" % fname)
   

    #Now get the components/death_times and save output
    components, death_times, comp_sizes = oldest_components(data.T, Filt, BarCodes_Rips0, num_comp_to_keep)
    fname = nme + '%s-oldest-comps.npy' % num_comp_to_keep 
    np.save(fname, components)
    fname = nme + '%s-oldest-comps_deathtimes.txt' % num_comp_to_keep 
    np.savetxt(fname, np.array(death_times))
    fname = nme + '%s-oldest-comps_sizes.npy' % num_comp_to_keep 
    np.save(fname, comp_sizes)
    print(' Components and deathtimes have been saved.')

    
    print(100*'-')
    print('       LONG-LIVED LOOPS:')

    #Save the loops birth/death times
    fname_births = nme + '%s-oldest-loops_birthtimes.txt' % num_comp_to_keep 
    fname_deaths = nme + '%s-oldest-loops_deathtimes.txt' % num_comp_to_keep
    np.savetxt(fname_births, np.array(loop_births))
    np.savetxt(fname_deaths, np.array(loop_deaths))
    print(' Birth/death times of loops have been saved.')

    
    #Write the loops to file
    print(' Producing persloop compatible input data...')
    fmn = nme+'f.txt'
    with open(fmn,'w') as f:
        f.write("%d %d\n"%(data.shape[1],data.shape[0]))
        #for data_slice in data:
        #    print(data_slice)
        np.savetxt(f, data, fmt='%-f',delimiter=' ')
            # Writing out a break to indicate different slices...
        #    f.write('\n')
        idx = data.shape[0]
        for j,ft in enumerate(Filt):
            if len(ft[0]) > 1 and len(ft[0]) < 4:
                #print(ft)
                f.write("# %d\n"%idx)
                f.write("i "+" ".join(str(i) for i in ft[0])+"\n")
                idx+=1
                Filt[j] = ft + (idx,)

    fmn = nme + "pers"

    #This section turns barcodes indexed by distance
    #into barcodes indexed by events in the filtration
    #
    bcode = []
    with open(fmn,'w') as f:
        for l in BarCodes_Rips0:
            if l[0] == 1:
                birth_time = -1
                if l[1][0] == 0.:
                    birth_time = 1
                else:
                    bfound = False
                    for ft in Filt:
                        if np.abs(l[1][0] - ft[1]) < peps:
                            birth_time = ft[2]
                            bfound = True
                            break
                    if not bfound:
                        print(" Warning, not found birth_time index!")
                death_time = np.inf
                if l[1][1] != np.inf:
                    dfound = False
                    for ft in Filt:
                        if np.abs(l[1][1] - ft[1]) < peps:
                            death_time = ft[2]
                            dfound = True
                            break
                    if not dfound:
                        print(" Warning, not found death_time index!")
                if death_time == np.inf:
                    death_time = 'inf'
                else:
                    death_time ='%i'%death_time
                f.write("%d %i %s\n"%(l[0],birth_time,death_time))
                bcode.append([l[0],birth_time,death_time])

    return Filt,bcode

def printloops(BarCodeObj):
    
    birth_times = []
    death_times = []
    for i in BarCodeObj:
        if i[0] == 1:
            birth = i[1][0]
            death = i[1][1]
            birth_times.append(birth)
            death_times.append(death)
            print(' Birth dist %5.3f, Death dist %5.3f'%(birth, death))

    return birth_times, death_times


def rmfloops(nme):
    import shutil
    from pathlib import Path
    import shutil
    dirpath = Path(nme+'floops')
    if dirpath.exists() and dirpath.is_dir():
        print(' Old floops directory found, removing...')
        try:
            shutil.rmtree(dirpath)
        except:
            print(" Warning, error when removing floops, check the permissions")
    return 

def runpersloop(nme):
    #Execute persloop, to be run after proc4persloop
    import subprocess
  
    print(' Starting persloop routine')
    rmfloops(nme)
    if threeD:
        exc = f"{persloopfolder}/persloop"
    else:
        #This should be the version of persloop compiled
        #from persloop-src-all-dim
        exc = f"{persloopfolder}/persloop-src-all-dim/build/persloop"
    command = '%s -f %sf.txt -s %spers'%(exc,nme,nme)
    print(' Executing: %s'%command)
    result = subprocess.check_output([command], stderr=subprocess.STDOUT,shell=True)
    output = result.decode("utf-8")
    print(' It is written')
    #    status = os.system(command)
    #if status == 0:
    #    print(' Persloop successfully ran')
    print("--------------------------------------")
    print("         Output:")
    print("--------------------------------------")
    print(output)
    if "off" in output[-4:]:
        print(" Persloop successfully ran")
    return 0
    

def getloop(num,dr):
    #Load a loop from a dir
    fn = dr+'/%i.off'%num
    data = []
    with open(fn,'r') as f:
        for i,l in enumerate(f.readlines()):
            if i == 1:
                tmp =l.split(" ")[0]
                nele = int(tmp)
                #nele = ([int(float(j)) for j in l.split(" ")][0])                
            if i > 1 and (i <=nele+1):# or allD): #Note this doesn't work for 6D
            #if i > 1: #For >3D
                tmp  =l.split(" ")
                if len(tmp) > 1:
                    data.append([float(j) for j in l.split(" ") if j not in ['\n'," "]])
    data = np.array(data)
    return data

def loopnums(dr):
    #Get all the loop files in dir
    import os
    import re
    contents = os.listdir(dr)# +'floops')
    loops = []
    for f in contents:
        procname = re.match('([0-9]+)\.off',f)
        if procname != None:
            if procname.group(0) == f:
                loops.append(int(procname.group(1)))
    return loops
