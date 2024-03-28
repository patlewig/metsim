import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib.request, urllib.parse, json
import time
from rdkit import Chem
import datetime
import requests
import pprint
import json
from metsim.utl.hcd import metsim_hcd_out


def metsim_run_toolbox_api(host_name = None, tb_port = None,
                       simulator_num = 15,
                       smiles = None,
                       casrn = None,
                       dtxsid = None,
                       chem_name = None,
                       idx = None):
    import datetime
    import time
    import pandas as pd
    import urllib.request, urllib.parse, json
    #simulator num is 15 for in vitro and 8 for in vivo
    metsim_url_base = f"""http://{host_name}:{tb_port}/api/v6/Metabolism/"""
    with urllib.request.urlopen(metsim_url_base) as url:
        oecd_metsim_guids = json.loads(url.read().decode())
    #store GUID of metabolism simulator
    guid = oecd_metsim_guids[simulator_num]['Guid']
    #Store base metsim info in output dictionary:
    oecd_dict = {'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                 'software': 'OECD QSAR Toolbox WebAPI',
                 'version': 6,
                 'params':{'depth': 3,
                           'organism': 'Rat',
                           'site_of_metabolism': False,
                           'model': [oecd_metsim_guids[simulator_num]['Caption']]
                          }
                }

    #make lambda function to perform metsim in API:
    metsim_exec = lambda url: json.loads(urllib.request.urlopen(url).read().decode())
    #Lambda function to search with url-encoded SMILES:
    search_base = f"""http://{host_name}:{tb_port}/api/v6/Search/"""
    search_smiles = lambda smiles: json.loads(urllib.request.urlopen(search_base+'smiles/false/true?smiles='+smiles).read().decode())
    #Lambda function to search on casrn if qsar_ready_smiles = nan:
    search_cas = lambda casrn: json.loads(urllib.request.urlopen(search_base+'cas/'+casrn+'/true').read().decode())
    stereo_filter = ['@','/','\\','.']
    oecd_dict['input'] = {'smiles': smiles,
                          'inchikey': None,
                          'casrn': casrn,
                          'hcd_smiles': None,
                          'dtxsid': dtxsid,
                          'chem_name': chem_name
                         }
    oecd_metab_list = []
    if pd.notna(smiles):
        #url encode SMILES
        smiles_encoded = urllib.parse.quote_plus(smiles)
        #Try metsim with base SMILES call first before doing anything more complicated than that:
        if pd.notna(smiles_encoded):
            print('Attempting metsim from SMILES input for index #'+str(idx)+'...')
            metsim_url = metsim_url_base+guid+'?smiles='+smiles_encoded  
            oecd_metab_list = metsim_exec(metsim_url) #store metabolite list for the input precursor
        if len(oecd_metab_list) > 0:
            oecd_dict['output'] = []
            oecd_dict['output'].append({'precursor': oecd_dict['input'],
                                        'successors': [{'enzyme': None,
                                                        'mechanism': None,
                                                        'metabolite': {'smiles': oecd_metab_list[j],
                                                                       'inchikey': None,
                                                                       'casrn': None,
                                                                       'hcd_smiles': None,
                                                                       'dtxsid': None,
                                                                       'chem_name': None
                                                                      },
                                                       } for j in range(len(oecd_metab_list))]
                                      })
            print('metsim succeeded for index #'+str(idx))
        else:
            oecd_dict = idx
    elif pd.notna(casrn) & ('NOCAS' not in casrn):
        oecd_dict = idx
    else:
        #This dictionary returns if no SMILES or CASRN are given:
        oecd_dict['input'] = {'smiles': None,
                              'inchikey': None,
                              'casrn': None,
                              'hcd_smiles': None,
                              'dtxsid': None,
                              'chem_name': None
                             }
        oecd_dict['output'] = [{'precursor': {'smiles': None,
                                              'inchikey': None,
                                              'casrn': None,
                                              'hcd_smiles': None,
                                              'dtxsid': None,
                                              'chem_name': None
                                             },
                                'successors': [{'enzyme': [],
                                                'mechanism': None,
                                                'metabolite': {'smiles': None,
                                                               'inchikey': None,
                                                               'casrn': None,
                                                               'hcd_smiles': None,
                                                               'dtxsid': None,
                                                               'chem_name': None
                                                               }
                                              }]
                              }]
        print('metsim failed for index #'+str(idx)+'. Neither SMILES nor CASRN were provided.')
    return (idx, oecd_dict)






    
def metsim_tb_search_logkow(casrn = None, host_name = None, tb_port = None, idx = None):
    """ 
    Search the OECD Toolbox database via its WebAPI for a chemical ID for an input chemical, and then 
    return the octanol-water partition coefficient to have a measure of its hydrophobicity.
    
    Inputs: 
    casrn: CAS Registry Number
    tb_port: Port number selected for locally running instance of the Toolbox Server
    
    Outputs:
    log_kow: Log10 scaled octanol-water partition coefficient, if available.
    """
    import pandas as pd
    import urllib.request, urllib.parse, json
    
    search_base = f"""http://{host_name}:{tb_port}/api/v6/search/"""
    kow_url = f"""http://{host_name}:{tb_port}/api/v6/calculation/41552380-4d5d-4eab-bee0-03774c0eabb6/"""
    search_cas = lambda casrn: json.loads(urllib.request.urlopen(search_base+'cas/'+casrn+'/true').read().decode())
    kow_exec = lambda chem_Id: json.loads(urllib.request.urlopen(kow_url+chem_Id).read().decode())
    stereo_filter = ['@','/','\\','.']
    chem_entries = []
    if casrn == None:
        print('No CASRN provided.')
        return None
    try:
        if 'NOCAS' not in casrn:
            print('Searching Toolbox database for ChemIds using CASRN: '+casrn+'.')
            cas_nohyphen = ''.join(casrn.split('-'))
            chem_entries = search_cas(cas_nohyphen)
        else:
            print('Invalid CASRN (contains "NOCAS"). Cannot search for ChemIds using this identifier.')
            return None
    except:
        print('Bad http request for index #'+str(idx)+' CASRN')
        return None
    if len(chem_entries) > 0:
        print('Toolbox database ChemIds found for CASRN: '+casrn+' via Toolbox API search.')
        #remove chem name lists from dicts so that duplicate dicts can be removed via set comprehension
        rem_names = [chem_entries[i].pop('Names',None) for i in range(len(chem_entries))]
        chem_entries = [dict(t) for t in {tuple(chem_entries[i].items()) for i in range(len(chem_entries))}] #remove duplicate entries
        #Filter out mixtures, and substances with no casrn:
        chem_mono = [chem_entries[j] for j in range(len(chem_entries)) if (chem_entries[j]['SubstanceType'] == 'MonoConstituent' and chem_entries[j]['Cas'] != 0)]
        #Store ChemID so long as smiles has no stereo information:
        chem_Id = [chem_mono[j]['ChemId'] for j in range(len(chem_mono)) if sum([stereo_filter[k] in chem_mono[j]['Smiles'] for k in range(len(stereo_filter))]) == 0]
        if len(chem_Id) > 0:
            print('Monoconstituent ChemIds found within search results...')
            for i in range(len(chem_Id)):
                print('Calculating Log Kow value for ChemId '+str(i+1)+'/'+str(len(chem_Id))+'...')
                try:
                    log_kow = kow_exec(chem_Id[i])
                    if log_kow != None:
                        if log_kow['Value'] != None:
                            print('Log Kow value successfully determined for CASRN: '+casrn+'.')
                            return float(log_kow['Value'])
                        else:
                            print('Log Kow value not found for current ChemId.')
                            continue
                except:
                    print('Bad http request for index #'+str(idx)+' Log Kow.')
                    return None
            if log_kow != None:
                if log_kow['Value'] == None:
                    print('ChemId(s) found, but no Log Kow available for CASRN: '+casrn+'.')
                    return None
            else:
                print('ChemId(s) found, but no Log Kow available for CASRN: '+casrn+'.')
                return None
        else:
            print('No ChemId(s) found for CASRN: '+casrn+'.')
            return None


def metsim_search_toolbox_api(tb_port = 16384,
                              host_name = None,
                              simulator_num = 15,
                              smiles = None,
                              casrn = None,
                              dtxsid = None,
                              chem_name = None,
                              idx = None):
    #Version of toolbox_metsim_api that does not do the SMILES metsim query, but searches for ChemIds on SMILES and CASRN to do metsim.
    #Run serially to circumvent TB search server crash issues until rectified by LMC updates.
    import datetime
    import time
    import pandas as pd
    import urllib.request, urllib.parse, json
    
    metsim_url_base = f"""http://{host_name}:{tb_port}/api/v6/Metabolism/"""
    with urllib.request.urlopen(metsim_url_base) as url:
        oecd_metsim_guids = json.loads(url.read().decode())
    #store GUID of metabolism simulator
    guid = oecd_metsim_guids[simulator_num]['Guid']
    #Store base metsim info in output dictionary:
    oecd_dict = {'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                 'software': 'OECD QSAR Toolbox WebAPI',
                 'version': 6,
                 'params':{'depth': 3,
                           'organism': 'Rat',
                           'site_of_metabolism': False,
                           'model': [oecd_metsim_guids[simulator_num]['Caption']]
                          }
                }

    #make lambda function to perform metsim in API:
    metsim_exec = lambda url: json.loads(urllib.request.urlopen(url).read().decode())
    #Lambda function to search with url-encoded SMILES:
    search_base = f"""http://{host_name}:{tb_port}/api/v6/search/"""
    search_smiles = lambda smiles: json.loads(urllib.request.urlopen(search_base+'smiles/false/true?smiles='+smiles).read().decode())
    #Lambda function to search on casrn if qsar_ready_smiles = nan:
    search_cas = lambda casrn: json.loads(urllib.request.urlopen(search_base+'cas/'+casrn+'/true').read().decode())
    stereo_filter = ['@','/','\\','.'] #comment if using daylight smiles with stereochemistry, uncomment for qsar-ready smiles
    oecd_dict['input'] = {'smiles': smiles,
                          'inchikey': None,
                          'casrn': casrn,
                          'hcd_smiles': None,
                          'dtxsid': dtxsid,
                          'chem_name': chem_name
                         }
    oecd_metab_list = []
    smiles_encoded = None
    if pd.notna(smiles):
        #url encode SMILES
        smiles_encoded = urllib.parse.quote_plus(smiles)
    if pd.notna(casrn):
        print('Base SMILES query for index #'+str(idx)+' yields no metabolites. Searching for alternate ChemIds...')
        #casrn without hyphens:
        if 'NOCAS' not in casrn:
            cas_nohyphen = ''.join(casrn.split('-'))
        chem_entries = []
        #Query database on both smiles and casrn within try statements.
        #Necessary becasue of HTTP Errors for some chemicals:
        if pd.notna(smiles_encoded):
            try: 
                chem_entries = search_smiles(smiles_encoded)
            except:
                print('Bad http request for index #'+str(idx)+' SMILES')
                chem_entries = []
        if len(chem_entries) > 0:
            try:
                if 'NOCAS' not in casrn:
                    chem_cas = search_cas(cas_nohyphen)
                    chem_entries = chem_entries+chem_cas
            except:
                print('Bad http request for index #'+str(idx)+' CASRN')
        else:
            try:
                if 'NOCAS' not in casrn:
                    chem_entries = search_cas(cas_nohyphen)
            except:
                print('Bad http request for index #'+str(idx)+' CASRN')
                chem_entries = []
        if len(chem_entries) > 0:
            #remove chem name lists from dicts so that duplicate dicts can be removed via set comprehension
            rem_names = [chem_entries[i].pop('Names',None) for i in range(len(chem_entries))]
            chem_entries = [dict(t) for t in {tuple(chem_entries[i].items()) for i in range(len(chem_entries))}] #remove duplicate entries
            #Filter out mixtures, and substances with no casrn:
            chem_mono = [chem_entries[j] for j in range(len(chem_entries)) if (chem_entries[j]['SubstanceType'] == 'MonoConstituent' and chem_entries[j]['Cas'] != 0)]
            #Store ChemID so long as smiles has no stereo information:
            chem_Id = [chem_mono[j]['ChemId'] for j in range(len(chem_mono)) if sum([stereo_filter[k] in chem_mono[j]['Smiles'] for k in range(len(stereo_filter))]) == 0] #uncomment if using qsar-ready smiles
            # chem_Id = [chem_mono[j]['ChemId'] for j in range(len(chem_mono))] #uncomment if using daylight smiles with stereochemistry
            if len(chem_Id) > 0:
                print('Alternate ChemId corresponding to a QSAR-Ready SMILES found for index #'+str(idx)+'. Attempting metsim...')
                #Cosntruct URL from base call for metsims + GUID of the simulator + ChemID:
                #In most cases, anything filtered down this far likely has the same SMILES 
                #even if ChemID is different, will yield same metabolite list.
                #So just take the first one.
                metsim_url = metsim_url_base+guid+'/'+chem_Id[0]
                oecd_metab_list = metsim_exec(metsim_url)
            elif len(chem_mono) > 0:
                print('Metsim failed. Alternate monoconstituent ChemIds found. Attempting metsim for alternate ChemId 1/'+str(len(chem_mono))+' for index #'+str(idx)+'...')
                #If the only monoconstituent db entries have a stereo SMILES associated with our search, 
                #just run it using that ChemId, will likely yield metabolites over an empty list.
                chem_Id = chem_mono[0]['ChemId']
                metsim_url = metsim_url_base+guid+'/'+chem_Id
                oecd_metab_list = metsim_exec(metsim_url)
            if len(oecd_metab_list) > 0:
                oecd_dict['output'] = []
                oecd_dict['output'].append({'precursor': oecd_dict['input'],
                                            'successors': [{'enzyme': None,
                                                            'mechanism': None,
                                                            'metabolite': {'smiles': oecd_metab_list[j],
                                                                           'inchikey': None,
                                                                           'casrn': None,
                                                                           'hcd_smiles': None,
                                                                           'dtxsid': None,
                                                                           'chem_name': None
                                                                          },
                                                           } for j in range(len(oecd_metab_list))]
                                          })
                print('Metsim succeeded for index #'+str(idx))
            elif len(chem_mono) > 1:
                print('Metsim failed. Attempting metsim from alternate monoconstituent ChemIds for index #'+str(idx)+'...')
                #If first chemID yields no metabolites, and there are other monoconstituent 
                #Chem_Ids to try run them until a non-empty metabolite list is yielded
                for l in range(1,len(chem_mono)):
                    oecd_metab_list = metsim_exec(metsim_url_base+guid+'/'+chem_mono[l]['ChemId'])
                    if len(oecd_metab_list) > 0:
                        oecd_dict['output'] = []
                        oecd_dict['output'].append({'precursor': oecd_dict['input'],
                                                    'successors': [{'enzyme': None,
                                                                    'mechanism': None,
                                                                    'metabolite': {'smiles': oecd_metab_list[j],
                                                                                   'inchikey': None,
                                                                                   'casrn': None,
                                                                                   'hcd_smiles': None,
                                                                                   'dtxsid': None,
                                                                                   'chem_name': None
                                                                                  },
                                                                   } for j in range(len(oecd_metab_list))]
                                                  })
                        print('Metsim succeeded for index #'+str(idx))
                        break
                    else:
                        #No available ChemIDs yield a metabolite list, inspect list manually later.
                        oecd_dict['output'] = [{'precursor': oecd_dict['input'],
                                                'successors': [{'enzyme': [],
                                                                'mechanism': None,
                                                                'metabolite': {'smiles': None,
                                                                               'inchikey': None,
                                                                               'casrn': None,
                                                                               'hcd_smiles': None,
                                                                               'dtxsid': None,
                                                                               'chem_name': None
                                                                              }
                                                              }]
                                              }]
                        print('Metsim failed for index #'+str(idx)+'. Alternate ChemId '+str(l+1)+'/'+str(len(chem_mono))+' yielded no metabolites.')
            else:
                oecd_dict['output'] = [{'precursor': oecd_dict['input'],
                                        'successors': [{'enzyme': [],
                                                        'mechanism': None,
                                                        'metabolite': {'smiles': None,
                                                                       'inchikey': None,
                                                                       'casrn': None,
                                                                       'hcd_smiles': None,
                                                                       'dtxsid': None,
                                                                       'chem_name': None
                                                                      }
                                                      }]
                                      }]
                print('Metsim failed for index #'+str(idx)+'. No alternate ChemIds found.')            
        else:
            oecd_dict['output'] = [{'precursor': oecd_dict['input'],
                                    'successors': [{'enzyme': [],
                                                    'mechanism': None,
                                                    'metabolite': {'smiles': None,
                                                                   'inchikey': None,
                                                                   'casrn': None,
                                                                   'hcd_smiles': None,
                                                                   'dtxsid': None,
                                                                   'chem_name': None
                                                                  }
                                                  }]
                                  }]
            print('metsim failed for index #'+str(idx)+'. Neither SMILES nor CASRN yield valid ChemIds.')
    return (idx, oecd_dict)


