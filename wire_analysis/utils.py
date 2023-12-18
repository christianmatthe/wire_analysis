# Collection of functions  that might be needed in other scripts
import json

def save_json_dict(dict_path, dict):
    with open((dict_path), 'w', encoding='utf-8') as f:
        json.dump(dict, f, ensure_ascii=False, indent=4)
    return

def load_json_dict(dict_path):
    with open(dict_path, 'r', encoding='utf-8') as f:
        dict_load = json.load(f)
    return dict_load