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
    Used to convert an input SMILES string into QSAR-Ready SMILES. Returns InChIKey structural identifier as well,
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
    Includes "hcd_smiles" as output dictionary key containing QSAR-Ready version of the input SMILES.
    
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
