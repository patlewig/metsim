
import os
import subprocess
import numpy as np
import pandas as pd
import urllib.request, urllib.parse, json
import time
import tempfile
from rdkit import Chem
import datetime
#import multiprocess as mp
from itertools import compress
import pprint
from metsim.utl.hcd import metsim_hcd_out



def metsim_metadata_full(metsim_out = [], fnam = None, metsim_cache = None):
    
    if len(metsim_out) > 0:
        if metsim_cache != None:
            #Supplement metadata via serial HCD query through individual input chemicals, precursors, successors/metabolites for a full metsim dataset
            for i in range(len(metsim_out)): # i = number of input chemicals
                if metsim_out[i]['input']['inchikey'] != None:
                        continue
                if metsim_out[i]['input']['smiles'] not in [cache_item['smiles'] for cache_item in metsim_cache]:
                    metsim_out[i]['input'] = metsim_hcd_out(smiles = metsim_out[i]['input']['smiles'],
                                                            casrn = metsim_out[i]['input']['casrn'],
                                                            dtxsid = metsim_out[i]['input']['dtxsid'],
                                                            chem_name = metsim_out[i]['input']['chem_name'])
                    metsim_cache.append(metsim_out[i]['input'])
                    print('Input query added to metadata cache...')
                else:
                    print('Input SMILES found in cached results. Inserting into dictionary...')
                    metsim_out[i]['input'] = metsim_cache[[idx for idx in range(len(metsim_cache)) if metsim_cache[idx]['smiles'] == metsim_out[i]['input']['smiles']][0]]
                for j in range(len(metsim_out[i]['output'])): # j = number of unique precursors
                    if 'likelihood' in list(metsim_out[i]['output'][j]['precursor'].keys()):
                        if metsim_out[i]['output'][j]['precursor']['smiles'] not in [cache_item['smiles'] for cache_item in metsim_cache]:
                            metsim_out[i]['output'][j]['precursor'] = metsim_hcd_out(smiles = metsim_out[i]['output'][j]['precursor']['smiles'],
                                                                                     casrn = metsim_out[i]['output'][j]['precursor']['casrn'],
                                                                                     dtxsid = metsim_out[i]['output'][j]['precursor']['dtxsid'],
                                                                                     chem_name = metsim_out[i]['output'][j]['precursor']['chem_name'],
                                                                                     likely = metsim_out[i]['output'][j]['precursor']['likelihood'])
                            metsim_cache.append(metsim_out[i]['output'][j]['precursor'])
                            print('Precursor query added to metadata cache...')
                        else:
                            print('Precursor SMILES found in cached results. Inserting into dictionary...')
                            metsim_out[i]['output'][j]['precursor'] = metsim_cache[[idx for idx in range(len(metsim_cache)) if metsim_cache[idx]['smiles'] == metsim_out[i]['output'][j]['precursor']['smiles']][0]]
                    else:
                        if metsim_out[i]['output'][j]['precursor']['smiles'] not in [cache_item['smiles'] for cache_item in metsim_cache]:
                            metsim_out[i]['output'][j]['precursor'] = metsim_hcd_out(smiles = metsim_out[i]['output'][j]['precursor']['smiles'],
                                                                                     casrn = metsim_out[i]['output'][j]['precursor']['casrn'],
                                                                                     dtxsid = metsim_out[i]['output'][j]['precursor']['dtxsid'],
                                                                                     chem_name = metsim_out[i]['output'][j]['precursor']['chem_name'])
                            metsim_cache.append(metsim_out[i]['output'][j]['precursor'])
                            print('Precursor query added to metadata cache...')
                        else:
                            print('Precursor SMILES found in cached results. Inserting into dictionary...')
                            metsim_out[i]['output'][j]['precursor'] = metsim_cache[[idx for idx in range(len(metsim_cache)) if metsim_cache[idx]['smiles'] == metsim_out[i]['output'][j]['precursor']['smiles']][0]]
                    for k in range(len(metsim_out[i]['output'][j]['successors'])): # k = number of metabolites per precursor
                        if 'likelihood' in list(metsim_out[i]['output'][j]['successors'][k]['metabolite'].keys()):
                            if metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles'] not in [cache_item['smiles'] for cache_item in metsim_cache]:
                                metsim_out[i]['output'][j]['successors'][k]['metabolite'] = metsim_hcd_out(smiles = metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles'],
                                                                                                           casrn = metsim_out[i]['output'][j]['successors'][k]['metabolite']['casrn'],
                                                                                                           dtxsid = metsim_out[i]['output'][j]['successors'][k]['metabolite']['dtxsid'],
                                                                                                           chem_name = metsim_out[i]['output'][j]['successors'][k]['metabolite']['chem_name'],
                                                                                                           likely = metsim_out[i]['output'][j]['successors'][k]['metabolite']['likelihood'])
                                metsim_cache.append(metsim_out[i]['output'][j]['successors'][k]['metabolite'])
                                print('Successor metabolite query added to metadata cache...')
                            else:
                                print('Successor metabolite SMILES found in cached results. Inserting into dictionary...')
                                metsim_out[i]['output'][j]['successors'][k]['metabolite'] = metsim_cache[[idx for idx in range(len(metsim_cache)) if metsim_cache[idx]['smiles'] == metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles']][0]] 
                        else:
                            if metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles'] not in [cache_item['smiles'] for cache_item in metsim_cache]:
                                metsim_out[i]['output'][j]['successors'][k]['metabolite'] = metsim_hcd_out(smiles = metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles'],
                                                                                                           casrn = metsim_out[i]['output'][j]['successors'][k]['metabolite']['casrn'],
                                                                                                           dtxsid = metsim_out[i]['output'][j]['successors'][k]['metabolite']['dtxsid'],
                                                                                                           chem_name = metsim_out[i]['output'][j]['successors'][k]['metabolite']['chem_name'])
                                metsim_cache.append(metsim_out[i]['output'][j]['successors'][k]['metabolite'])
                                print('Successor metabolite query added to metadata cache...')
                            else:
                                print('Successor metabolite SMILES found in cached results. Inserting into dictionary...')
                                metsim_out[i]['output'][j]['successors'][k]['metabolite'] = metsim_cache[[idx for idx in range(len(metsim_cache)) if metsim_cache[idx]['smiles'] == metsim_out[i]['output'][j]['successors'][k]['metabolite']['smiles']][0]] 
                        print('input: '+str(i+1)+'/'+str(len(metsim_out))+' precursor: '+str(j+1)+'/'+str(len(metsim_out[i]['output']))+' metabolite: '+str(k+1)+'/'+str(len(metsim_out[i]['output'][j]['successors'])))
                if fnam != None:
                    json.dump(metsim_out, open(fnam,'w'))
        else:
            return metsim_metadata_full(metsim_out = metsim_out, fnam = fnam, metsim_cache = [])
    else:
        raise('Please supply a metsim dataset (list of dictionaries)')
    # print(metsim_out)
    return metsim_out    
    
#for generational tracking of output:
def recursive_gen_list(input_df = None,
                       parent_list = [],
                       successor_list = None,
                       out_list = [],
                       gen = 1):
    if type(input_df) == 'NoneType':
        raise('Please include BioTransformer output csv dataframe as input_df.')
    successor_dict = {}
    for i in parent_list.index:
        if pd.isna(parent_list[i]):
            successor_df = input_df.loc[pd.isna(input_df['Precursor ID']),:]
            successor_list = successor_df['Metabolite ID']
            print(str(len(successor_list))+' Successors found for parent chemical')
        else:
            successor_df = input_df.loc[input_df['Precursor ID'] == parent_list[i],:]
            successor_list = successor_df['Metabolite ID']
        if len(successor_list) > 0:

            successor_dict = {'precursor': {'smiles': successor_df['Precursor SMILES'].unique()[0],
                                        'inchikey': None,
                                        'casrn': None,
                                        'hcd_smiles': None,
                                        'dtxsid': None,
                                        'chem_name': None
                                       },
                          'successors': [{'enzyme': successor_df.loc[j,'Enzyme(s)'].split('\n'),
                                          'mechanism': successor_df.loc[j,'Reaction'],
                                          'generation': gen,
                                          'metabolite': {'smiles': successor_df.loc[j,'SMILES'],
                                                         'inchikey': None,
                                                         'casrn': None,
                                                         'hcd_smiles': None,
                                                         'dtxsid': None,
                                                         'chem_name': None
                                                        }
                                       } for j in successor_df.index]
                             }
            if pd.notna(parent_list[i]):
                print(str(len(successor_df))+' successors found for gen '+str(gen)+' precursor '+parent_list[i])
            if successor_dict['precursor']['smiles'] in [out_list[j]['precursor']['smiles'] for j in range(len(out_list))]:
                print('Precursor-successor relationship already recorded, skipping to next precursor.')
            else:
                out_list = out_list+[successor_dict]
                print('dict added to out_list for index '+str(i))
                #check for more generations:
                for k in successor_list.index:
                    print('starting recursion on successor index '+str(k))
                    out_list = recursive_gen_list(input_df = input_df,
                                                  parent_list = pd.Series(successor_list[k]),
                                                  successor_list = [],
                                                  out_list = out_list,
                                                  gen = gen+1)
                    print('recursion complete for successor index '+str(k))
            return out_list
        else:
            print('No successors for gen '+str(gen)+' precursor '+parent_list[i]+', moving to next gen 1 precursor.')
            gen = 1
            return out_list


# Function 3: BioTransformer 3.0 Human Phase I + Phase II Metsim.
def btrans_metsim_subprocess(models = ['cyp450','phaseII'],
                             cyp_mode = 1,
                             cycles = [2,1],
                             smiles = None,
                             casrn = None,
                             del_tmp = False,
                             dtxsid = None,
                             chem_name = None,
                             idx = None,
                             multi_proc = False):
    
    '''
    
    input: SMILES, List of Models, List of Cycles, Delete Tempfile (True/False)
action: Simulates human metabolism (Default: 2 cycles Phase I/"cyp450", 1 cycle Phase II/"phaseII") using BioTransformer 3.0 through Java in the command prompt, producing output data in temporary files that can be kept or deleted. Recursively searches through output CSV to process data into a standardized dictionary output.
output: Tuple with dummy index (for parallel processing), Dictionary of precursor and successor SMILES, CASRN, DTXSID, InChIKey as supplemented by HCD (or RDKit for InChIKey), and filename (if del_tmp = False).

    '''
    
# Need to determine how to implement recursion with multiprocessing...
    
    #set up initial dictionary outputs:
    btrans_dict = {'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                   'software': 'BioTransformer',
                   'version': 3.0,
                   'params':{'depth': sum([cycles[i] if ('allHuman' != models[i]) and ('superbio' != models[i]) else 4 for i in range(len(cycles))]),
                             'organism': 'Human',
                             'site_of_metabolism': False,
                             'model': [str(cycles[i])+'x '+models[i] for i in range(len(models))]
                            }
                  }
    btrans_dict['input'] = {'smiles': smiles,
                            'inchikey': None,
                            'casrn': casrn,
                            'hcd_smiles': None,
                            'dtxsid': dtxsid,
                            'chem_name': chem_name
                           } 
    
    #change directory to BioTransformer directory
    btrans_dir = 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/GenRA/BioTransformer3.0'
    # btrans_dir = 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/GenRA/BioTransformer1.1.5'
    if os.curdir != btrans_dir:
        os.chdir(btrans_dir)
       
    if pd.notna(smiles): #check that we have a SMILES string input
        #Begin constructing BioTransformer input command string:
        btrans_cmd = 'java -jar BioTransformer3.0_20220615.jar -k pred ' #current version
        # btrans_cmd = 'java -jar BioTransformer-1.1.5.jar -k pred ' #2019 paper version
        # btrans_cmd = 'java -jar BioTransformer-1-0-6.jar -k pred ' #original version (Only uses single models!)
        #incorporate model parameters into command string:
        if len(models) == 1:
            btrans_model = '-b "'+models[0]+'" '
            if cycles[0] > 1:
                btrans_model = btrans_model+'-s '+str(cycles[0])
        else:
            btrans_model = '-q "'
            for i in range(len(models)):
                if i < len(models)-1:
                    btrans_model = btrans_model+models[i]+':'+str(cycles[i])+';'
                else:
                    btrans_model = btrans_model+models[i]+':'+str(cycles[i])+'" '
        btrans_cmd = btrans_cmd+btrans_model+' -cm '+str(cyp_mode)+' ' #specify cyp_mode with later versions of BioTransformer
        # btrans_cmd = btrans_cmd+btrans_model+' ' #no cyp_mode in input command for original 2019 relase of BioTransformer
        
        #create temporary file to store BioTransformer output:
        print('Generating output tempfile for index '+str(idx)+', DTXSID: '+str(dtxsid)+'...')
        fnum, fnam = tempfile.mkstemp(dir = './tmpfiles',
                                      prefix = 'btrans_out_'+dtxsid+'_',
                                      suffix = '.csv')
        
        #Finish constructing input command using the tempfile name:
        fnam_split = fnam.split('\\')
        btrans_fnam = '.\\'+fnam_split[-2]+'\\'+fnam_split[-1]
        #insert input smiles string and output csv tempfile name:
        btrans_cmd = btrans_cmd + '-ismi "'+smiles+'" -ocsv "'+btrans_fnam+'"' #use Daylight SMILES input
        # btrans_cmd = btrans_cmd + '-isdf "'+smiles+'" -ocsv "'+btrans_fnam+'"'
        print(btrans_cmd)
        #Run BioTransformer with subprocess and the appropriate input command string:
        print('beginning metsim for index #'+str(idx)+'...')
        btrans_out = subprocess.run(btrans_cmd,
                                    shell=True,
                                    stdout=subprocess.PIPE) 
        
        #Preallocate data frame and read tempfile output if BioTransformer ran successfully:
        btrans_out_df = pd.DataFrame()
        if btrans_out.stderr == None:
            try:
                btrans_out_df = pd.read_csv(fnam)
            except:
                pass
        else:
            print('error, check stdout')
            print(btrans_out.stdout.decode())
            
        btrans_dict['output'] = []
        if multi_proc == False:
            #Process metabolism data (if it exists for the given smiles):
            if os.path.getsize(fnam) > 0:
                input_df = pd.read_csv(fnam)
                parent_list = input_df['Precursor ID']
                parent_list = parent_list.drop_duplicates()
                btrans_dict['output'] = recursive_gen_list(input_df = input_df,
                                                           parent_list = parent_list,
                                                           successor_list = [],
                                                           out_list = [],
                                                           gen = 1)
            else: #Valid SMILES given, no metabolites produced:
                print('No metabolites produced for index #'+str(i))
                btrans_dict['output'] = [{'precursor': btrans_dict['input'],
                                          'successors': [{'enzyme': [],
                                                          'mechanism': None,
                                                          'generation': None,
                                                          'metabolite': {'smiles': None,
                                                                         'inchikey': None,
                                                                         'casrn': None,
                                                                         'hcd_smiles': None,
                                                                         'dtxsid': None,
                                                                         'chem_name': None
                                                                        }
                                                        }]
                                        }]
        else:
            #Need to determine how to multiprocess recursion, until then, store input parameters and output structure.
            btrans_dict['output'] = [{'precursor': btrans_dict['input'],
                                      'successors': [{'enzyme': [],
                                                      'mechanism': None,
                                                      'generation': None,
                                                      'metabolite': {'smiles': None,
                                                                     'inchikey': None,
                                                                     'casrn': None,
                                                                     'hcd_smiles': None,
                                                                     'dtxsid': None,
                                                                     'chem_name': None
                                                                    }
                                                    }]
                                    }]
        #close temporary output file. Will delete if del_tempfile function input is set to True.
        os.close(fnum)
        if del_tmp == True:
            os.remove(fnam)
    else:
        print('No SMILES string provided for index #'+str(idx)+".")
        #This list returns if no metabolites are formed/BioTransformer fails.
        btrans_dict['output'] = [{'precursor': {'smiles': None,
                                                'inchikey': None,
                                                'casrn': None,
                                                'hcd_smiles': None,
                                                'dtxsid': None,
                                                'chem_name': None
                                               },
                                  'successors': [{'enzyme': [],
                                                  'mechanism': None,
                                                  'generation': None,
                                                  'metabolite': {'smiles': None,
                                                                 'inchikey': None,
                                                                 'casrn': None,
                                                                 'hcd_smiles': None,
                                                                 'dtxsid': None,
                                                                 'chem_name': None
                                                                }
                                                }]
                                }]
        btrans_dict['input'] = {'smiles': None,
                                'inchikey': None,
                                'casrn': None,
                                'hcd_smiles': None,
                                'dtxsid': None,
                                'chem_name': None
                               }
    print('MetSim complete for index '+str(idx)+'.')
    if del_tmp == False:
        #until figuring out how to implement recursive metabolite search with multiprocessing, return filename for separate analysis function.
        return (idx, btrans_dict, fnam)
    else:
        return (idx, btrans_dict, None)    
    
    
    

