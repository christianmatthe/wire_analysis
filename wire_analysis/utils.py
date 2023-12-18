# Collection of functions  that might be needed in other scripts
import json
import numpy as np

def save_json_dict(dict_path, dict):
    with open((dict_path), 'w', encoding='utf-8') as f:
        json.dump(dict, f, ensure_ascii=False, indent=4)
    return

def load_json_dict(dict_path):
    with open(dict_path, 'r', encoding='utf-8') as f:
        dict_load = json.load(f)
    return dict_load

def load_extractor_dict_json(
                                dict_path
                                ):
    extractor_dict_unsorted_lists = load_json_dict(dict_path=dict_path)
    extractor_dict_unsorted  ={}
    for key,val in extractor_dict_unsorted_lists.items():
        # every numpy array must become a list
        extractor_dict_unsorted[key] = (np.array(
            extractor_dict_unsorted_lists[key]))
    return extractor_dict_unsorted