import os
import numpy as np
import pandas as pd
import time
import tempfile
from rdkit import Chem
import datetime
from metsim.utl.hcd import metsim_hcd_out


def metsim_run_times(times_df = None,
                     model = 'vitro'):
    
    """
    Takes the metabolism report loaded as a pandas dataframe from the CSV output from the TIMES Desktop Application for a given simulator and returns the metsim hierarchy
    with all available metadata from EPA Cheminformatics Modules API for all precursor-successor relationships predicted from either the In Vivo Rat Simulator 
    or the In Vitro Rat Liver S9 model within times. 
    
    How to get the CSV required:
    In the tab named after the simulator in the TIMES GUI, generate the metabolism "Map" from the "Reports" toolbar after running the simulator on all parent chemicals,
    The "Map Report" table that is produced can be copied to the clipboard, pasted into spreadsheet software, and saved as a CSV file. Load this as a dataframe via Pandas.
    
    Inputs:
    times_df (dataframe, required): Dataframe of metabolism output from TIMES
    model (str, optional): Specify either In Vitro Rat Liver S9 with "vitro" (default), or In Vivo Rat Simulator with "vivo".
    
    Output:
    times_metsim_dict (list of dict): Generationally tracked, MetSim hierarchy structured dictionary of precursor-successor relationships for all available input chemicals 
    """
    
    
    
    
    if times_df.empty == False:
        precursor_df = times_df.loc[pd.notna(times_df['Chemical name']),:]
        precursor_df = precursor_df.drop_duplicates(subset='Chemical name', ignore_index=True)
        successor_df = times_df.loc[pd.isna(times_df['Chemical name']),:]
        #print(successor_df)
        times_metsim_dict = []
        times_cache = []
        if model.lower() == 'vivo':
            model = 'In Vivo Rat Simulator v08.14'
        if model.lower() == 'vitro':
            model = 'In Vitro Rat Liver S9 v12.18'
        for i in range(len(precursor_df)):
            # print('i = ',i)
            if precursor_df.loc[precursor_df.index[i],'SMILES'] not in [cache_item['smiles'] for cache_item in times_cache]:
                times_metsim_dict.append({'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                                          'software': 'OASIS TIMES',
                                          'version': '2.31.2.82',
                                          'params':{'depth': 3,
                                                    'organism': 'Rat',
                                                    'site_of_metabolism': False,
                                                    'model': [model]
                                                   },
                                          'input': metsim_hcd_out(smiles = precursor_df.loc[precursor_df.index[i],'SMILES'], dtxsid = precursor_df.loc[precursor_df.index[i],'Chemical name']),
                                          'output': []
                                        })
                times_cache.append(times_metsim_dict[-1]['input'])
                print('Input query added to metadata cache...')
            else:
                print('Input SMILES found in cached results. Inserting into dictionary...')
                times_metsim_dict.append({'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                                          'software': 'OASIS TIMES',
                                          'version': '2.31.2.82',
                                          'params':{'depth': 3,
                                                    'organism': 'Rat',
                                                    'site_of_metabolism': False,
                                                    'model': [model]
                                                   },
                                          'input': times_cache[[idx for idx in range(len(times_cache)) if times_cache[idx]['smiles'] == precursor_df.loc[precursor_df.index[i],'SMILES']][0]],
                                          'output': []
                                        })
            tmpdf = successor_df.loc[successor_df['Parent Number'] == precursor_df.loc[precursor_df.index[i],'Parent Number'],:]
            for j in range(len(tmpdf['Predecessor ID'].unique())):
                # print('j = ',j)
                times_metsim_dict[i]['output'].append({'precursor': None, 
                                                       'successors': []
                                                     })
                for k  in range(len(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], :])):
                    # print('k = ',k)
                    k_idx = tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], :].index
                    if tmpdf.loc[tmpdf.index[j], 'Predecessor ID'] != str(1):
                        if tmpdf.loc[tmpdf['ID of metabolite'] == int(tmpdf.loc[tmpdf.index[j],'Predecessor ID']), 'SMILES'].values[0] not in [cache_item['smiles'] for cache_item in times_cache]:
                            times_metsim_dict[i]['output'][j]['precursor'] = metsim_hcd_out(tmpdf.loc[tmpdf['ID of metabolite'] == int(tmpdf.loc[tmpdf.index[j],'Predecessor ID']), 'SMILES'].values[0])
                            times_cache.append(times_metsim_dict[-1]['output'][-1]['precursor'])
                            print('precursor metabolite query added to metadata cache...')
                        else:
                            print('Precursor SMILES found in cached results. Inserting into dictionary...')
                            times_metsim_dict[i]['output'][j]['precursor'] = times_cache[[idx for idx in range(len(times_cache)) if times_cache[idx]['smiles'] == tmpdf.loc[tmpdf['ID of metabolite'] == int(tmpdf.loc[tmpdf.index[j],'Predecessor ID']), 'SMILES'].values[0]][0]]
                        
                        if tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]] not in [cache_item['smiles'] for cache_item in times_cache]:
                            times_metsim_dict[i]['output'][j]['successors'].append({'enzyme': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation type'][k_idx[k]],
                                                                                    'mechanism': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation name'][k_idx[k]],
                                                                                    'metabolite': metsim_hcd_out(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]]),
                                                                                    'generation': int(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Level of generation'][k_idx[k]])
                                                                                   })
                            times_cache.append(times_metsim_dict[-1]['output'][-1]['successors'][-1]['metabolite'])
                            print('successor metabolite query added to metadata cache...')
                        else:
                            print('Successor SMILES found in cached results. Inserting into dictionary...')
                            times_metsim_dict[i]['output'][j]['successors'].append({'enzyme': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation type'][k_idx[k]],
                                                                                    'mechanism': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation name'][k_idx[k]],
                                                                                    'metabolite': times_cache[[idx for idx in range(len(times_cache)) if times_cache[idx]['smiles'] == tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]]][0]],
                                                                                    'generation': int(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Level of generation'][k_idx[k]])
                                                                                   })
                    else:
                        times_metsim_dict[i]['output'][j]['precursor'] = times_metsim_dict[i]['input']
                        if tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]] not in [cache_item['smiles'] for cache_item in times_cache]:
                            times_metsim_dict[i]['output'][j]['successors'].append({'enzyme': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation type'][k_idx[k]],
                                                                                    'mechanism': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation name'][k_idx[k]],
                                                                                    'metabolite': metsim_hcd_out(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]]),
                                                                                    'generation': int(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Level of generation'][k_idx[k]])
                                                                                   })
                            times_cache.append(times_metsim_dict[-1]['output'][-1]['successors'][-1]['metabolite'])
                            print('successor metabolite query added to metadata cache...')
                        else:
                            print('Successor SMILES found in cached results. Inserting into dictionary...')
                            times_metsim_dict[i]['output'][j]['successors'].append({'enzyme': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation type'][k_idx[k]],
                                                                                    'mechanism': tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Transformation name'][k_idx[k]],
                                                                                    'metabolite': times_cache[[idx for idx in range(len(times_cache)) if times_cache[idx]['smiles'] == tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'SMILES'][k_idx[k]]][0]],
                                                                                    'generation': int(tmpdf.loc[tmpdf['Predecessor ID'] == tmpdf.loc[tmpdf.index[j],'Predecessor ID'], 'Level of generation'][k_idx[k]])
                                                                                   })
    return times_metsim_dict
