'''
Combine mesa outputs (history + frequencies), and store in hdf5 format.
Three versions:
1)  Complete grid: 2+3.
2)  Classic grid: parameters that are present in the ".history" files.
    From protostar to the end (Dnu~10muHz),
    including every time step.
3)  Radial frequency grid: parameters in the ".history" files, 
    as well as radial frequencies in the "sum" files.
    From the MS turnoff to the end (Dnu~10muHz),
    including once every 5 time step.
'''

import numpy as np
import os
import sys
from astropy.io import ascii
from astropy.table import Table, Column
import corner
import h5py
import re
import pandas as pd




rootpath = '/project/RDS-FSC-astero-RW/low-mass-red-giants/'
ifmulti, Nthread = True, int(sys.argv[1])

sys.path.append(rootpath)

from lib.toolkit import history, sums



if __name__ == '__main__':

    work_dir = rootpath + 'hpc/'

    coarse_grid = pd.read_csv(work_dir+'sapphire/template/sapphire_grid_input_params.txt', sep=',\s+', engine='python', index_col='index')
    # coarse_grid2 = pd.read_csv(work_dir+'coarse_v1/template/coarse_grid_input_params_v1_add.txt', sep=',\s+', engine='python', index_col='index')
    # coarse_grid = pd.concat((coarse_grid1,coarse_grid2))
    coarse_grid['index'] = coarse_grid.index

    # retrive a list for all history files
    historyPaths = np.array([work_dir+'sapphire/heb/history/'+f for f in os.listdir(work_dir+'sapphire/heb/history/') if f.endswith('.history')])
    h1 = np.array([f for f in os.listdir(work_dir+'sapphire/heb/history/') if f.endswith('.history')])
    h2 = np.array([f.split('.h5')[0] for f in os.listdir(work_dir+'sapphire/heb/complete_grid/') if f.endswith('.h5')])
    idx = np.isin(h1, h2)
    historyPaths = historyPaths[~idx]

    # sumPaths = [] 
    # subfolders = [work_dir+'coarse_v1/freqs/'+f+'/' for f in os.listdir(work_dir+'coarse_v1/freqs/') if f.startswith('index')]
    # for subfolder in subfolders:
    #     sumPaths += [subfolder+f for f in os.listdir(subfolder) if f.endswith('.FGONG.sum')]
    # sumDirs = np.array([f.split('/index')[0]+'/' for f in sumPaths])
    # sumNames = np.array([f.split('/')[-1] for f in sumPaths])


    # for historyPath in historyPaths[:]:
    def multi(historyPath):

        historyFile = re.split('/',historyPath)[-1]
        historyIndexStr = re.split('.history', re.split('index',historyFile)[-1])[0]
        historyIndex = int(historyIndexStr)

        outdirCompleteGrid = work_dir+'sapphire/heb/complete_grid/'

        # # read in models
        profileIndexPath = historyPath.split('.history')[0]+'profile.index'
        if os.path.exists(profileIndexPath):
            h = history(historyPath, ifReadProfileIndex=True)
            ifprofile = True 
        else:
            h = history(historyPath, ifReadProfileIndex=False)
            ifprofile = False 

        track = h.track 
        _, idx = np.unique(track['model_number'], return_index=True)
        track = track[idx]

        table = pd.DataFrame(track)

        if ifprofile:
            table = table.merge(pd.DataFrame(h.profileIndex), on='model_number', how='left')

        # # grid initial parameters
        idx = coarse_grid['index']==historyIndex
        index, Xinit, Yinit, Zinit, amlt, fov_shell, fov0_shell, fov_core, fov0_core = coarse_grid.iloc[np.where(idx)[0][0]][['index', 'Xinit', 'Yinit', 'Zinit', 'amlt', 'fov_shell', 'fov0_shell', 'fov_core', 'fov0_core']].to_list()
        table['index'], table['Xinit'], table['Yinit'], table['Zinit'], table['amlt'] = index, Xinit, Yinit, Zinit, amlt
        table['fov_shell'], table['fov0_shell'], table['fov_core'], table['fov0_core'] = fov_shell, fov0_shell, fov_core, fov0_core


        # # assign a phase
        phase = np.zeros(len(table))
        turnoffidx = table['center_h1'] < 1.e-7
        msidx = (10.0**(table['log_Lnuc']-table['log_L'])>0.99) & (~turnoffidx)
        pmsidx = (10.0**(table['log_Lnuc']-table['log_L'])<=0.99) & (~turnoffidx)
        sgidx = (turnoffidx) & (table['nu_max']>=300)
        rgbidx = (turnoffidx) & (table['nu_max']<300)
        hebidx = (turnoffidx) & ((table['center_he4']+table['center_he3']) <0.95)
        tfidx = (turnoffidx) & (~sgidx) & (~rgbidx) & (~hebidx)

        phase[pmsidx] = -1
        phase[msidx] = 0
        phase[sgidx] = 1
        phase[rgbidx] = 2
        phase[hebidx] = 3
        phase[tfidx] = -2
        table['phase'] = phase

        # # log properties
        table['luminosity'] = 10.0**table['log_L']
        table['radius'] = 10.0**table['log_R']
        table['Teff'] = 10.0**table['log_Teff']

        # # seismic scaling quantities
        Dnu_sun, numax_sun, Teff_sun = 135.1, 3090., 5777.
        table['delta_nu_scaling'] = table['star_mass']**0.5 * table['radius']**-1.5 * Dnu_sun
        table['numax_scaling'] = table['star_mass'] * table['radius']**-2.0 * (table['Teff']/Teff_sun)**-0.5 * numax_sun
        
        # # surface quantities
        Zsun, Xsun = 0.0134, 0.7381 # 0.0134, 0.7381, a09 # 0.0169, 0.7345, gs98
        table['FeH'] = np.log10((1-table['surface_h1']-table['surface_he4']-table['surface_he3'])/table['surface_h1']) - np.log10(Zsun/Xsun)


        # # # assign a prior
        # age = np.array(table['star_age'])
        # prior = np.concatenate([[0],(age[2:]-age[:-2])/2.,[0]])
        # prior[~np.isfinite(prior)] = 0. 
        # prior = prior/np.sum(prior)
        # table['prior'] = prior


        subfolder = work_dir+'sapphire/heb/freqs/index{:06.0f}/'.format(historyIndex)
        sumPaths = [subfolder+f for f in os.listdir(subfolder) if f.endswith('.FGONG.sum')]
        sumDirs = np.array([f.split('/index')[0]+'/' for f in sumPaths])
        sumNames = np.array([f.split('/')[-1] for f in sumPaths])

        if ifprofile:
            # # set a seismic flag
            table['flag_seismo'] = np.array(np.isfinite(table['profile_number']), dtype=int)

            # # read in radial mode frequencies
            seismicCols = ['l', 'n_p', 'n_g', 'n_pg', 'E', 'E_p', 'E_g', 'E_norm', 'freq']
            seismicData = [[] for i in range(len(seismicCols))]

            for imod, mod in table.loc[:,:].iterrows():
                profileIndex = mod['profile_number']
                if (not np.isfinite(profileIndex)):
                    for i in range(len(seismicData)):
                        seismicData[i].append(np.nan)
                else:
                    sumFile = 'index{:06.0f}profile{:0.0f}.data.FGONG.sum'.format(historyIndex, profileIndex)
                    if len(sumDirs[sumNames==sumFile])==0:
                        for i in range(len(seismicData)):
                            seismicData[i].append(np.nan)
                        table.loc[imod, 'flag_seismo'] = 0
                        table.loc[imod, 'profile_number'] = np.nan
                    else:
                        sumPath = sumDirs[sumNames==sumFile][0] + 'index{:06.0f}/'.format(historyIndex) + sumFile
                        s = sums(sumPath, verbose=False)
                        for i in range(len(seismicData)-1):
                            seismicData[i].append(np.array([s.modeSummary[seismicCols[i]]]))
                        seismicData[-1].append(np.array([s.modeSummary['Refreq']]))
        else:
            # # set a seismic flag
            table['flag_seismo'] = np.zeros(len(table), dtype=int)

        # #  write out the table
        with h5py.File(outdirCompleteGrid+'{:s}.h5'.format(historyFile), 'w') as h5f:
            for col in table.columns:
                h5f.create_dataset(col, data=table[col])
            if ifprofile:
                for imod, mod in table.iterrows():
                    profileIndex = mod['profile_number']
                    if (not np.isfinite(profileIndex)): continue
                    for i in range(len(seismicCols)):
                        h5f.create_dataset('profile{:0.0f}/{:s}'.format(profileIndex, seismicCols[i]), data=seismicData[i][imod])

        
        return 0.



    if ifmulti:
        from multiprocessing import Pool 
        with Pool(Nthread) as p:
            res = p.map(multi, historyPaths)
    else:
        for historyPath in historyPaths:
            multi(historyPath)


