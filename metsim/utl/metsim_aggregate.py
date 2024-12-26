import datetime

def metsim_aggregate_results(metsim_list = None, metsim_aggregate_out = None):
    """
    Function to aggregate two or more MetSim hierarchies together into one aggregate hierarchy that specifies how many models predicted each metabolite,
    and which models predicted each metabolite.
    
    Args:
        metsim_list (required): processed MetSim hierarchy. Must contain at least SMILES and inchikey metadata. 
                                *note: If unprocessed predictions from a tool get inserted, will raise errors. To process data, 
                                use metsim_metadata_full with your unprocessed predictions to generate the proper input for metsim_list.

        metsim_aggregate_out (optional): If an aggregate set of predictions already exists, and adding another metsim_list to it is desired, 
                                         input the existing processed aggregate hierarchy of predictions.
    """
    if metsim_list != None:
        
        for i in range(len(metsim_list)):
            #if the aggregate list is empty, insert the first tool's output for the first chemical and update the output to reflect which simulator was used, 
            #and that--so far, one tool predicted these metabolites:
            if metsim_aggregate_out == None:
                metsim_aggregate_out = []
                metsim_aggregate_out.append({'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                                            'software': 'Ensemble Metabolic Predictions',
                                            'input': metsim_list[i]['input'],
                                            'output': metsim_list[i]['output']})
                for j in range(len(metsim_aggregate_out[0]['output'])):
                    for k in range(len(metsim_aggregate_out[0]['output'][j]['successors'])):
                        if 'toolbox' in metsim_list[i]['software'].lower():
                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['toolbox_vivo_rat_phase1']
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['toolbox_vitro_rat_phase1']
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'times' in metsim_list[i]['software'].lower():
                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['times_vivo_rat_phase1_phase2']
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['times_vitro_rat_phase1_phase2']
                                metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'biotransformer' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['biotransformer_human_phase1_phase2']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'admet' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['admet_human_phase1']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'epa' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['cts_human_phase1']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'acd' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['predicted_by'] = ['acdlabs_percepta']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'smpdb' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['observed_in'] = ['smpdb']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_sources_observed'] = 1
                        if 'jcim' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['observed_in'] = ['jcim_literature']
                            metsim_aggregate_out[0]['output'][j]['successors'][k]['num_sources_observed'] = 1
            #logical test to see if the ith parent chemical data already exists in the aggregate dictionary.
            #If not present, insert it into the list:
            elif sum([metsim_list[i]['input']['chem_name'] in metsim_aggregate_out[j]['input']['chem_name'] for j in range(len(metsim_aggregate_out))]) == 0:
                metsim_aggregate_out.append({'datetime': str(datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')),
                                            'software': 'Aggregate Metabolic Predictions',
                                            'input': metsim_list[i]['input'],
                                            'output': metsim_list[i]['output']})
                print('lenth of pred dict aggregate ',len(metsim_aggregate_out))
                for j in range(len(metsim_aggregate_out[-1]['output'])):
                    for k in range(len(metsim_aggregate_out[-1]['output'][j]['successors'])):
                        if 'toolbox' in metsim_list[i]['software'].lower():
                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['toolbox_vivo_rat_phase1']
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['toolbox_vitro_rat_phase1']
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'times' in metsim_list[i]['software'].lower():
                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['times_vivo_rat_phase1_phase2']
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['times_vitro_rat_phase1_phase2']
                                metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'biotransformer' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['biotransformer_human_phase1_phase2']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'admet' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['admet_human_phase1']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'epa' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['cts_human_phase1']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'acd' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['predicted_by'] = ['acdlabs_percepta']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_models_predicted'] = 1
                        if 'smpdb' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['observed_in'] = ['smpdb']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_sources_observed'] = 1
                        if 'jcim' in metsim_list[i]['software'].lower():
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['observed_in'] = ['jcim_literature']
                            metsim_aggregate_out[-1]['output'][j]['successors'][k]['num_sources_observed'] = 1
            #All parent chemicals accounted for by at least one tool. Check to see if the tool from "metsim_list" is represented in the output of the aggregate pathway. If not, add those outputs:
            #Check precursors first. If precursor isn't present in the existing precursors for that parent, add it and its successors.
            #If precursors are all accounted for, add successors for each tool, append "predicted_by" list if another tool also predicts those successors, and +1 "num_tools_predicted"
            else:
                if sum([metsim_list[i]['input']['chem_name'] in metsim_aggregate_out[j]['input']['chem_name'] for j in range(len(metsim_aggregate_out))]) != 0:
                    idx_in = [j for j in range(len(metsim_aggregate_out)) if metsim_list[i]['input']['chem_name'] == metsim_aggregate_out[j]['input']['chem_name']][0]
                    input_precursors_inchi2d = [metsim_list[i]['output'][j]['precursor']['inchikey'].split('-')[0] for j in range(len(metsim_list[i]['output']))]
                    aggregate_precursors_inchi2d = [metsim_aggregate_out[idx_in]['output'][j]['precursor']['inchikey'].split('-')[0] for j in range(len(metsim_aggregate_out[idx_in]['output']))]
                    for j in range(len(input_precursors_inchi2d)):
                        #For the ith parent, if the jth output precursor from a tool's predictions isn't present in the existing precursors for the parent, append the jth precursor and its successors.
                        if sum([input_precursors_inchi2d[j] in aggregate_precursors_inchi2d[k] for k in range(len(aggregate_precursors_inchi2d))]) == 0:
                            metsim_aggregate_out[idx_in]['output'].append(metsim_list[i]['output'][j])
                            #Iterate over successors for the newly appended output to reflect the tool predicted by, and number of tools predicted.
                            for k in range(len(metsim_aggregate_out[idx_in]['output'][-1]['successors'])):
                                if 'toolbox' in metsim_list[i]['software'].lower():
                                    if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['toolbox_vivo_rat_phase1']
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                    if 's9' in metsim_list[i]['params']['model'][0].lower():
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['toolbox_vitro_rat_phase1']
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'times' in metsim_list[i]['software'].lower():
                                    if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['times_vivo_rat_phase1_phase2']
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                    if 's9' in metsim_list[i]['params']['model'][0].lower():
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['times_vitro_rat_phase1_phase2']
                                        metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'biotransformer' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['biotransformer_human_phase1_phase2']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'admet' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['admet_human_phase1']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'epa' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['cts_human_phase1']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'acd' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['predicted_by'] = ['acdlabs_percepta']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_models_predicted'] = 1
                                if 'smpdb' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['observed_in'] = ['smpdb']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_sources_observed'] = 1
                                if 'jcim' in metsim_list[i]['software'].lower():
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['observed_in'] = ['jcim_literature']
                                    metsim_aggregate_out[idx_in]['output'][-1]['successors'][k]['num_sources_observed'] = 1
                        #Else, for the ith parent, if the jth output precursor from a tool's predictions already exists in the list of precursors,
                        #Check that the output from that tool isn't already inventoried in the aggregate
                        #If not, incorporate its successors to that output, or add to the number of tools predicted and predicted by lists.
                        elif sum([input_precursors_inchi2d[j] in aggregate_precursors_inchi2d[k] for k in range(len(aggregate_precursors_inchi2d))]) > 0:
                            #Passing this logical test says that the jth precursor is already in the aggregate list for the ith parent chemical.
                            #In case the length of the output lists is different due to appending new precursors into the aggregate list, find the index corresponding to the jth precursor for the tool in the aggregate precursor list:
                            idx_con = [k for k in range(len(metsim_aggregate_out[idx_in]['output'])) if metsim_aggregate_out[idx_in]['output'][k]['precursor']['inchikey'].split('-')[0] == metsim_list[i]['output'][j]['precursor']['inchikey'].split('-')[0]][0]
                            for k in range(len(metsim_list[i]['output'][j]['successors'])):
                                if metsim_list[i]['output'][j]['successors'][k]['metabolite']['inchikey'] != None:
                                    #If the kth successor of the jth precursor is not present in the aggregate list for that same precursor, add it into the successor list:
                                    if sum([metsim_list[i]['output'][j]['successors'][k]['metabolite']['inchikey'].split('-')[0] in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][l]['metabolite']['inchikey'].split('-')[0] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'])) if metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][l]['metabolite']['inchikey'] != None ]) == 0:
                                        metsim_aggregate_out[idx_in]['output'][idx_con]['successors'].append(metsim_list[i]['output'][j]['successors'][k])
                                        print('successor appended for ',j,'th precursor')
                                        if 'toolbox' in metsim_list[i]['software'].lower():
                                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['toolbox_vivo_rat_phase1']
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['toolbox_vitro_rat_phase1']
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'times' in metsim_list[i]['software'].lower():
                                            if 'vivo' in metsim_list[i]['params']['model'][0].lower():
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['times_vivo_rat_phase1_phase2']
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                            if 's9' in metsim_list[i]['params']['model'][0].lower():
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['times_vitro_rat_phase1_phase2']
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'biotransformer' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['biotransformer_human_phase1_phase2']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'admet' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['admet_human_phase1']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'epa' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['cts_human_phase1']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'acd' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['predicted_by'] = ['acdlabs_percepta']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_models_predicted'] = 1
                                        if 'smpdb' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['observed_in'] = ['smpdb']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_sources_observed'] = 1
                                        if 'jcim' in metsim_list[i]['software'].lower():
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['observed_in'] = ['jcim_literature']
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][-1]['num_sources_observed'] = 1
                                    else:
                                        #Successor is present in the successor list. Find its index in aggregate successors. Check which model corresponds to this prediction. 
                                        #if it's different from those inventoried currently, add it to the predicted by list, and increment number of models predicted
                                        metsim_aggregate_out[idx_in]['output'][idx_con]['successors'] = [metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'])) if metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][l]['metabolite']['inchikey'] != None]
                                        idx_suc_con = [l for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'])) if metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][l]['metabolite']['inchikey'].split('-')[0] == metsim_list[i]['output'][j]['successors'][k]['metabolite']['inchikey'].split('-')[0]][0]
                                        if 'toolbox' in metsim_list[i]['software'].lower():
                                            if 'vivo' in metsim_list[i]['params']['model'][0].lower() and sum(['vivo' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('toolbox_vivo_rat_phase1')
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                            if 's9' in metsim_list[i]['params']['model'][0].lower() and sum(['vitro' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('toolbox_vitro_rat_phase1')
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'times' in metsim_list[i]['software'].lower():
                                            if 'vivo' in metsim_list[i]['params']['model'][0].lower() and sum(['vivo' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('times_vivo_rat_phase1_phase2')
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                            if 's9' in metsim_list[i]['params']['model'][0].lower() and sum(['vitro' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('times_vitro_rat_phase1_phase2')
                                                metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'biotransformer' in metsim_list[i]['software'].lower() and sum(['biotransformer' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('biotransformer_human_phase1_phase2')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'admet' in metsim_list[i]['software'].lower() and sum(['admet' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('admet_human_phase1')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'epa' in metsim_list[i]['software'].lower() and sum(['epa' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0:
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('cts_human_phase1')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'acd' in metsim_list[i]['software'].lower() and sum(['acd' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by']))]) == 0: 
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['predicted_by'].append('acdlabs_percepta')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_models_predicted'] += 1
                                        if 'smpdb' in metsim_list[i]['software'].lower() and sum(['smpdb' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in']))]) == 0:
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in'].append('smpdb')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_sources_observed'] += 1
                                        if 'jcim' in metsim_list[i]['software'].lower() and sum(['jcim' in metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in'][l] for l in range(len(metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in']))]) == 0:
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['observed_in'].append('jcim_literature')
                                            metsim_aggregate_out[idx_in]['output'][idx_con]['successors'][idx_suc_con]['num_sources_observed'] += 1
    else:
        raise Exception("Please provide a list of MetSim predictions.")
    return metsim_aggregate_out