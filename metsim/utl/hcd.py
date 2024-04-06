import urllib.request
import os
import pandas as pd
import numpy as np
import requests
import time
from rdkit import Chem
import datetime
import pprint
import json


def metsim_hcd_out(smiles = None, 
                   casrn = None,
                   dtxsid = None,
                   chem_name = None,
                   likely = None):
    """
    Query function for the Cheminformatics Modules Standardizer API, formerly wrapped within the Hazard Comparison Dashboard (HCD) API. 
    Used to convert an input SMILES string into QSAR-Ready SMILES (hcd_smiles). Returns InChIKey structural identifier as well,
    along with any other chemical identifer metadata if available, and not already given as inputs (e.g., CASRN, DTXSID, Chemical Name).
    
    If SMILES is not known, but DTXSID is known, can instead query on DTXSID to obtain Daylight SMILES from the Comptox Chemicals Dashboard API (CCD API),
    and subsequently query the Standardizer API using the SMILES obtained from the CCD API.
    
    Required Inputs:
    smiles: Daylight SMILES string
    or
    dtxsid: DSSTox Substance Identifier
    
    Optional Inputs:
    chem_name: Chemical name, whether trade name or IUPAC
    casrn: Chemical Abstracts Services Registry Number
    inchikey: International Chemical Identifier Key (InChIKey)
    likely: If MetSim predictions are obtained from the Chemical Transformation Simulator, can optionally keep the transformation "likelihood" parameter
    
    Returns:
    out_dict: Output dictionary containing all available output data for the given chemical, using the input parameter names as dictionary keys.
    Includes "hcd" as output dictionary key containing QSAR-Ready version of the input SMILES.
    
    Examples:
    
    SMILES given as sole input:
    input:    
    test_dict = metsim_hcd_out(smiles = "OCCOCCO")
    
    output:
    Attempting query of Cheminformatics Modules Standardizer with SMILES: OCCOCCO...
    Query succeeded.
    test_dict
    {'smiles': 'OCCOCCO',
     'casrn': '111-46-6',
     'hcd_smiles': 'OCCOCCO',
     'inchikey': 'MTHSVFCYNBDYFN-UHFFFAOYNA-N',
     'dtxsid': 'DTXSID8020462',
     'chem_name': 'Diethylene glycol',
     'likelihood': None}
    
    DTXSID given as sole input:
    
    input:    
    test_dict = metsim_hcd_out(dtxsid = "DTXSID4020402")
    
    output:
    Attempting query of Comptox Chemicals Dashboard with DTXSID: DTXSID4020402...
    Query succeeded.
    No SMILES given. Using CCD output SMILES.
    Attempting query of Cheminformatics Modules Standardizer with SMILES: CC1=C(N)C=C(N)C=C1...
    Query succeeded.
    test_dict
    {'smiles': 'CC1=C(N)C=C(N)C=C1',
     'casrn': '95-80-7',
     'hcd_smiles': 'CC1C=CC(N)=CC=1N',
     'inchikey': 'VOZKAJLKRJDJLL-UHFFFAOYNA-N',
     'dtxsid': 'DTXSID4020402',
     'chem_name': '2,4-Diaminotoluene',
     'likelihood': None}
     
     Empty inputs:
     input:
     test_dict = metsim_hcd_out(smiles = None, dtxsid = None)
     
     output:
     test_dict
     {'smiles': None,
      'casrn': None,
      'hcd_smiles': None,
      'inchikey': None,
      'dtxsid': None,
      'chem_name': None,
      'likelihood': None}   
    """
    
    ccd_out = []
    if dtxsid != None and smiles == None:
        #get metadata from Comptox Chemicals Dashboard for a given DTXSID (No structure searching atm).
        ccd_url = 'https://comptox.epa.gov/dashboard-api/ccdapp2/chemical-detail/search/by-dsstoxsid?id='+dtxsid
        ccd_success = 0
        try_count = 0
        while ccd_success == 0 and try_count < 3:
            try:
                print('Attempting query of Comptox Chemicals Dashboard with DTXSID: '+dtxsid+'...')
                ccd_out = json.loads(urllib.request.urlopen(ccd_url).read().decode())
                
                ccd_success = 1
                try_count+=1
            except:
                #Given that this occasionally fails randomly due to timeout errors, 
                #but then works again later, try again after a 1 second pause.
                #Should work on second attempt.
                print('URL Error Occurred, reattempting CCD query in 0.5 seconds.')
                time.sleep(0.5)
                try_count+=1
            print('Query succeeded.')
    if smiles != None or len(ccd_out) > 0:
        if smiles != None:
            smiles_url = urllib.parse.quote_plus(smiles) #URL encode SMILES string.
        elif len(ccd_out) > 0 and ccd_out['smiles'] != None:
            print('No SMILES given. Using CCD output SMILES.')
            smiles = ccd_out['smiles']
            smiles_url = urllib.parse.quote_plus(smiles) #URL enconde CCD smiles string.
        else:
            print('No SMILES given, and no SMILES available from CCD output.')
            smiles_url = None
        
        base_url = "https://hcd.rtpnc.epa.gov/api/stdizer?workflow=qsar-ready&smiles=" #Production environment (current, no VPN needed)
        
        if smiles_url != None:
            hcd_url = base_url+smiles_url
            hcd_success = 0
            try_count = 0
            hcd_out = []
            while hcd_success == 0 and try_count < 3:
                try:
                    print('Attempting query of Cheminformatics Modules Standardizer with SMILES: '+smiles+'...')
                    time.sleep(0.5)
                    hcd_out = json.loads(urllib.request.urlopen(hcd_url).read().decode())
                    print('Query succeeded.')
                    hcd_success = 1
                    try_count+=1
                except:
                    #Given that this occasionally fails randomly due to timeout errors, 
                    #but then works again later, try again after a 0.5 second pause.
                    #Should work on second attempt.
                    print('URL Error Occurred, reattempting Cheminformatics Modules query in 0.5 seconds.')
                    time.sleep(0.5)
                    try_count+=1
            if len(hcd_out) > 0:
                out_dict = {'smiles': smiles, 
                            'casrn': casrn,
                            'hcd_smiles': hcd_out[0]['smiles'],
                            'inchikey': hcd_out[0]['inchiKey'],
                            'dtxsid': dtxsid,
                            'chem_name': chem_name,
                            'likelihood': likely}
                if out_dict['dtxsid'] == None:
                    if 'DTXSID' in hcd_out[0]['id']:
                        out_dict['dtxsid'] = hcd_out[0]['id']
                    elif len(ccd_out) > 0 and ccd_out['dsstoxSubstanceId'] != None:
                        out_dict['dtxsid'] = ccd_out['dsstoxSubstanceId']
                if out_dict['casrn'] == None:
                    if 'casrn' in hcd_out[0].keys():
                        out_dict['casrn'] = hcd_out[0]['casrn']
                    elif len(ccd_out) > 0 and ccd_out['casrn'] != None:
                        out_dict['casrn'] = ccd_out['casrn']
                if out_dict['chem_name'] == None:
                    if len(ccd_out) > 0 and ccd_out['preferredName'] != None:
                        out_dict['chem_name'] = ccd_out['preferredName']
                    elif 'name' in hcd_out[0].keys():
                        out_dict['chem_name'] = hcd_out[0]['name']
                if out_dict['inchikey'] == None and len(ccd_out) > 0:
                    if ccd_out['inchiKey'] != None:
                        out_dict['inchikey'] = ccd_out['inchiKey'] 
            else:
                out_dict = {'smiles': smiles,
                            'casrn': casrn,
                            'hcd_smiles': None,
                            'inchikey': None,
                            'dtxsid': dtxsid,
                            'chem_name': chem_name,
                            'likelihood': likely
                           }
                #HCD Returns empty list. Try to supplement with metadata from RDKit.
                try:
                    smiles_mol = Chem.MolFromSmiles(smiles)
                    out_dict['inchikey'] = Chem.inchi.MolToInchiKey(smiles_mol)
                    print('RDKit generated InChIKey for SMILES: ',smiles)
                except:
                    #Rarely, BioTransformer makes a bad SMILES string for a metabolite, and RDKit can't convert it to an InChIKey. Store None
                    print('RDKit failed to generate an inchikey for SMILES: '+smiles)
                    out_dict['inchikey'] = None
        else:
            out_dict = {'smiles': smiles,
                        'casrn': casrn,
                        'hcd_smiles': None,
                        'inchikey': None,
                        'dtxsid': dtxsid,
                        'chem_name': chem_name,
                        'likelihood': likely
                       }
    else:
        out_dict = {'smiles': smiles,
                    'casrn': casrn,
                    'hcd_smiles': None,
                    'inchikey': None,
                    'dtxsid': dtxsid,
                    'chem_name': chem_name,
                    'likelihood': likely
                   }
    return out_dict

def metsim_metadata_full(metsim_out = [], fnam = None, metsim_cache = None):
    """
    MetSim simulators return a hierarchically structured dictionary of predicted precursor and successor SMILES from the given input SMILES. This function
    Takes a completed dataset of MetSim predictions in list format, iteratively queries EPA Cheminformatics Modules Standardizer API via the single metsim_hcd_out
    utility function call for all precursors and successors in the metsim hierarchy, and returns the list with all available metadata for all the predicted precursors and metabolites.
    (i.e., each parent input smiles produces its own dictionary, start from the list of those dictionaries for a full dataset of input parent chemicals).
    
    incorporates temporary caching of the metsim_hcd_out queries to avoid repepetive queries for commonly reoccuring metabolites to save time.
    
    Input:
    
    metsim_out (list): list of MetSim hierarchically structured dictionary outputs for each input chemical in the dataset.
    fnam (str, optional): the filename you would like to save the output to, ending in .json (e.g., "my_output_data.json"), will save to current working directory if not otherwise specified in fnam.
    metsim_cache (list, optional): list of cached metsim_hcd_out queries.
    """
    if len(metsim_out) > 0:
        if type(metsim_out) == dict:
            metsim_out = [metsim_out] #if a single chemical dataset is given, convert to list of dict
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
                    json.dump(metsim_out, open(str(datetime.datetime.now().strftime('%Y-%m-%d'))+'_metsim_metadata_full.json','w'))
        else:
            return metsim_metadata_full(metsim_out = metsim_out, fnam = fnam, metsim_cache = [])
    else:
        raise('Please supply a metsim dataset (list of dictionaries)')
    # print(metsim_out)
    return metsim_out    
