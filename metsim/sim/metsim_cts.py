import pandas as pd
import numpy as np
import requests
import time
from rdkit import Chem
import datetime
import pprint
import json
from metsim.utl.hcd import metsim_hcd_out

def metsim_run_cts(smiles = None, gens = 0):
    """
    Using CTS ChemAxon Human Phase I metabolism via its REST API
    takes an input SMILES string that will serve as the value for the "structure" key in the POST JSON dict and an integer value for number of metabolism cycles (0-4)
    """
    
    if smiles != None:
        if gens < 5 and gens > 0:
            url = 'https://qed.epa.gov/cts/rest/metabolizer/run'
            post =  {"structure": smiles, "generationLimit": gens, "transformationLibraries": ["mammalian_metabolism"]}
            time.sleep(0.5)
            
            try:    
                output = requests.post(url = url, json = post, verify=False).json()
                print('CTS API successfully returned metabolism predictions for '+str(gens)+' cycles of Phase I metabolism for input SMILES: '+str(smiles))
                return output
            except:
                print('CTS API failed to generate metabolites.')
                return None
        else:
            print('Generations must be between 1 and 4.')
            return None
    else:
        return None
    
def metsim_cts_phaseI_out(in_smiles = None, cts_data = None, depth = 0, out_dict = None):
    """
    Run CTS Human Phase I metabolism simulator, then recursively process the CTS output into metsim hierarchy structure.
    inputs:
    in_smiles: input SMILES string
    depth: Transformation depth in terms of integer number of generations (1-4)
    cts_data: CTS output JSON dictionary
    out_dict: MetSim output dictionary
    """
    
    #generate initial output dictionary hierarchy, and obtain CTS JSON output from API:
    if out_dict == None:
        out_dict = {'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                     'software': 'EPA Chemical Transformation Simulator',
                     'version': '1.3.2.2',
                     'params':{'depth': depth,
                               'organism': 'human',
                               'site_of_metabolism': False,
                               'model': 'ChemAxon Human Phase I'
                              },
                     'input': {'smiles': in_smiles,
                               'inchikey': None,
                               'casrn': None,
                               'hcd_smiles': None,
                               'dtxsid': None,
                               'chem_name': None
                              },
                     'output': []
                   }
    
    #if no CTS output data exists, query the API
    if cts_data == None:            
        cts_data = metsim_run_cts(smiles = in_smiles, gens = depth)
        
        #if the API call fails and yields no output
        if cts_data == None or 'error' in list(cts_data.keys()):
            print('CTS API call yielded no output data.')
            out_dict['api_success'] = False
            successor_dict = {'precursor': {'smiles': out_dict['input']['smiles'],
                                            'inchikey': None,
                                            'casrn': None,
                                            'hcd_smiles': None,
                                            'dtxsid': None,
                                            'chem_name': None,
                                            'likelihood': None
                                           },
                               'successors': [{'enzyme': None,
                                               'mechanism': None,
                                               'generation': None,
                                               'metabolite': {'smiles': None,
                                                              'inchikey': None,
                                                              'casrn': None,
                                                              'hcd_smiles': None,
                                                              'dtxsid': None,
                                                              'chem_name': None,
                                                              'likelihood': None
                                                             }
                                              }]
                             }
            
            #At least supplement metadata with inchikey from RDKit if possible
            try:
                out_dict['input']['inchikey'] = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(out_dict['input']['smiles']))
                print('RDKit InChIKey generation for input SMILES: '+out_dict['input']['smiles']+' succeeded.')
                out_dict['output'].append(successor_dict) 
                return out_dict 
            except:
                print('RDKit InChIKey generation failed.')
                out_dict['output'].append(successor_dict) 
                return out_dict
            
    #successful API call will have a key in the output dictionary for "children"
    if 'children' in list(cts_data.keys()):
        print('children key found in dictionary')
        #if metabolites are present in this level, store their information
        if len(cts_data['children']) > 0:
            print(str(len(cts_data['children']))+' children found for precursor generation level '+str(cts_data['data']['generation']))
            out_dict['api_success'] = True
            successor_dict = {'precursor': {'smiles': cts_data['data']['smiles'],
                                            'inchikey': None,
                                            'casrn': None,
                                            'hcd_smiles': None,
                                            'dtxsid': None,
                                            'chem_name': None,
                                            'likelihood': cts_data['data']['likelihood']
                                           },
                               'successors': [{'enzyme': None,
                                               'mechanism': cts_data['children'][j]['data']['routes'],
                                               'generation': cts_data['children'][j]['data']['generation'],
                                               'metabolite': {'smiles': cts_data['children'][j]['data']['smiles'],
                                                              'inchikey': None,
                                                              'casrn': None,
                                                              'hcd_smiles': None,
                                                              'dtxsid': None,
                                                              'chem_name': None,
                                                              'likelihood': cts_data['children'][j]['data']['likelihood']
                                                             }
                                              } for j in range(len(cts_data['children']))]
                             }
            out_dict['output'] = out_dict['output']+[successor_dict]
            print('precursor-successor relationships appended for current generational level.')
            for i in range(len(cts_data['children'])):
                if len(cts_data['children'][i]['children']) > 0:
                    print('children found in next generational level, recursing...')
                    out_dict = metsim_cts_phaseI_out(in_smiles = None, cts_data = cts_data['children'][i], depth = depth, out_dict = out_dict)
            return out_dict
        else:
            #no children returned, but successful API call:
            print('CTS API call successful. No metabolites produced for SMILES: '+out_dict['input']['smiles'])
            out_dict['api_success'] = True
            successor_dict = {'precursor': {'smiles': cts_data['data']['smiles'],
                                            'inchikey': None,
                                            'casrn': None,
                                            'hcd_smiles': None,
                                            'dtxsid': None,
                                            'chem_name': None,
                                            'likelihood': None
                                           },
                               'successors': [{'enzyme': None,
                                               'mechanism': None,
                                               'generation': None,
                                               'metabolite': {'smiles': None,
                                                              'inchikey': None,
                                                              'casrn': None,
                                                              'hcd_smiles': None,
                                                              'dtxsid': None,
                                                              'chem_name': None,
                                                              'likelihood': None
                                                             }
                                              }]
                             }
            try:
                out_dict['output'][0]['precursor']['inchikey'] = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(out_dict['input']['smiles']))
                print('RDKit InChIKey generation for input SMILES: '+out_dict['input']['smiles']+' succeeded.')
                out_dict['output'].append(successor_dict) 
                return out_dict 
            except:
                print('RDKit InChIKey generation failed.')
                out_dict['output'].append(successor_dict) 
                return out_dict 
            return out_dict
    elif 'children' not in list(cts_data.keys()):
        cts_data = cts_data['data']
        return metsim_cts_phaseI_out(in_smiles = None, cts_data = cts_data, depth = depth, out_dict = out_dict)
