import os
import pandas
from os.path import join as pjoin
pfams_path = "proteom/pfams/"

def process_hmm_file(f) :
    domtblout_head = ["target_name" , "target_accession" , "query_name" , "query_accession" , "E-value","score-sequence" , "bias-sequence" , "bdE-value","score-best-domain" , "bias--best-domain" , "exp" , "reg" , "clu" , "ov" , "env" , "dom" , "rep" , "inc" , "description_of_target"]
    data = pandas.read_csv(pjoin(pfams_path,f), delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,18))
    data = data.loc[data['bdE-value'] < 10e-6]
    pfams_dict = {p : [] for p in set(data['target_name'])}
    for a,b in data.iterrows():
        pfams_dict[b['target_name']] += [b['query_name']]
    return pfams_dict
