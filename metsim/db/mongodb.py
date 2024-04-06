import pymongo
import pandas as pd
import json
import pprint
import os
from rdkit import Chem

def metsim_mongo_search(mongo_db = None, collection = None, inchikey = None, casrn = None, dtxsid = None, precursor = False):
    """
    PyMongo connection and search function for collections located within "metsim_v1" database on pb.epa.gov remote server.
    Searches specified collection within "metsim_v1" for records matching the given input chemical or structural identifiers.
    Search hierarchy prioritizes inchkey first, if inchikey is None, then casrn, and if casrn is None, then dtxsid.
    
    Returns full record if the record contains the desired identifier, specified in either the precursor/parent or successor/metabolite dictionary levels
    
    input parameters:
    -----------------
    mongo_db: the pymongo.database.Database object pointing to the mongo database of interest, (e.g., mongo_db = client['metsim_v1']).
    
    collection: string, name of the collection you wish to search in metsim_v1
    
    inchikey: string, InChIKey of the chemical you wish to search for. Partial matching available on 2D structure from first block of InChIKey.
    
    casrn: string, CAS Registry number of the chemical you wish to search for.
    
    dtxsid: string, DSSTox Structural ID for the chemical you wish to search for.
    
    precursor: boolean, specifying whether you wish to pull records containing precursors with the identifier of choice, or metabolites with identifier of choice.
    
    
    returns:
    --------
    list of dictionary records found by the PyMongo find() function. 
    
    Examples:
    Given initial connection to mongodb via pymongo:
    client =  pymongo.MongoClient(f"mongodb://{username}:{password}@pb.epa.gov/metsim_v1")
    db = client['metsim_v1']
    col = db['validation_jcheminfmodel_drugs_nsaids_phase1'] #J. Chem. Inf. Model. human metabolite data.
    insert data into collection via col.insert_one({metsim_dict}) or col.insert_many([metsim_dict1, metsim_dict2,...]) 
    
    test search cases:
    metsim_mongo_search(mongo_db = db, collection = col, precursor = False, inchikey = 'PGRMEHUQIPZHKM-UHFFFAOYNA-N') #find metabolite on inchikey
    metsim_mongo_search(mongo_db = db, collection = col, precursor = False, inchikey = 'PGRMEHUQIPZHKM-UHFFFAOYNA-P') #wrong charge to force inchikey first block search
    metsim_mongo_search(mongo_db = db, collection = col, precursor = False, casrn = '152405-02-2') #find metabolite on casrn
    metsim_mongo_search(mongo_db = db, collection = col, precursor = False, dtxsid = 'DTXSID00934492') #find metabolite on dtxsid 
    metsim_mongo_search(mongo_db = db, collection = col, precursor = True, inchikey = 'GIIZNNXWQWCKIB-UHFFFAOYNA-N') #find precursor on inchikey
    metsim_mongo_search(mongo_db = db, collection = col, precursor = True, inchikey = 'GIIZNNXWQWCKIB-UHFFFAOYNA-P') #wrong charge to force inchikey first block search
    metsim_mongo_search(mongo_db = db, collection = col, precursor = True, casrn = '89365-50-4') #find precursor on casrn
    metsim_mongo_search(mongo_db = db, collection = col, precursor = True, dtxsid = 'DTXSID6023571') #find precursor on dtxsid
    """
    #connect to metsim_v1 on pb.epa.gov and point to desired collection
    if type(mongo_db) != pymongo.database.Database:
        raise ValueError('"db" parameter must be of type pymongo.database.Database. Please use client = MongoClient("mongodb://<username>:<password>@pb.epa.gov/metsim_v1?authSource=admin") to connect to pb.epa.gov and specify db = client["metsim_v1"] and try again.')
    if collection in mongo_db.list_collection_names():
        mongo_col = mongo_db[collection]
        print('Accessing '+collection+' collection on pb.epa.gov/metsim_v1.')
    else:
        raise ValueError('"'+collection+'" collection does not exist in this database. Try a different collection name.')
    
    #check that precursor parameter is boolean
    if type(precursor) != type(True):
        raise ValueError('"precursor" can either be True or False.')
    
    #For metabolite search on either inchikey, casrn, or dtxsid
    if precursor == False:
        if pd.notna(inchikey):
            print('Searching '+collection+' for metabolites containing inchikey: '+inchikey+'...')
            search_list = list(mongo_col.find({'output.successors.metabolite.inchikey': inchikey}))
            if len(search_list) > 0:
                print('Found metabolite containing inchikey: '+inchikey+'.')
            else:
                print('No metabolites found for full inchikey. Searching again on 2D structural block: '+inchikey.split('-')[0]+'...')
                search_list = list(mongo_col.find({'output.successors.metabolite.inchikey': {'$regex':inchikey.split('-')[0]}}))
                if len(search_list) > 0:
                    print('Found metabolite containing inchikey 2D structural block: '+inchikey.split('-')[0]+'.')
                else:
                    print('No metabolite found containing inchikey: '+inchikey+'.')
            return search_list
        elif pd.notna(casrn):
            print('Searching '+collection+' for metabolites containing casrn: '+casrn+'...')
            search_list = list(mongo_col.find({'output.successors.metabolite.casrn': casrn}))
            if len(search_list) > 0:
                print('Found metabolite containing casrn: '+casrn+'.')
            else:
                print('No metabolite found containing casrn:'+casrn+'.')
            return search_list
        elif pd.notna(dtxsid):
            print('Searching '+collection+' for metabolites containing dtxsid: '+dtxsid+'...')
            search_list = list(mongo_col.find({'output.successors.metabolite.dtxsid': dtxsid}))
            if len(search_list) > 0:
                print('Found metabolite containing dtxsid: '+dtxsid+'.')
            else:
                print('No metabolite found containing dtxsid: '+dtxsid+'.')
            return search_list
        else:
            raise ValueError('Please give either an "inchikey", "casrn", or "dtxsid" as a search parameter.')
            
    #for precursor chemical search on either inchikey, casrn, or dtxsid        
    if precursor == True:
        if pd.notna(inchikey):
            print('Searching '+collection+' for precursor chemicals containing inchikey: '+inchikey+'...')
            search_list = list(mongo_col.find({'output.precursor.inchikey': inchikey}))
            if len(search_list) > 0:
                print('Found precursor chemical for inchikey: '+inchikey+'.')
            else:
                print('No precursor chemical found for full inchikey. Searching again on 2D structural block: '+inchikey.split('-')[0]+'...')
                search_list = list(mongo_col.find({'output.precursor.inchikey': {'$regex':inchikey.split('-')[0]}}))
                if len(search_list) > 0:
                    print('Found precursor chemical containing inchikey 2D structural block: '+inchikey.split('-')[0]+'.')
                else:
                    print('No precursor chemical found containing inchikey: '+inchikey+'.')
            return search_list
        elif pd.notna(casrn):
            print('Searching '+collection+' for precursor chemical from casrn input...')
            search_list = list(mongo_col.find({'output.precursor.casrn': casrn}))
            if len(search_list) > 0:
                print('Found precursor chemical for the given casrn.')
            else:
                print('No precursor chemical found for the given casrn.')
            return search_list
        elif pd.notna(dtxsid):
            print('Searching '+collection+' for precursor chemials from dtxsid input...')
            search_list = list(mongo_col.find({'output.precursor.dtxsid': dtxsid}))
            if len(search_list) > 0:
                print('Found precursor chemical containing dtxsid: '+dtxsid+'.')
            else:
                print('No precursor chemical found containing dtxsid: '+dtxsid+'.')
            return search_list
        else:
            raise ValueError('Please give either an "inchikey", "casrn", or "dtxsid" as a search parameter.')
