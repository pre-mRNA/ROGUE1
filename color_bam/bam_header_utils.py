import os
import subprocess

# fetch ROGUE1 version 
def get_version():
    try:
        cwd = os.path.dirname(os.path.abspath(__file__))
        repo_root = subprocess.check_output(['git', '-C', cwd, 'rev-parse', '--show-toplevel']).decode('utf-8').strip()
        
        # fetch latest git tag
        git_tag = subprocess.check_output(['git', '-C', repo_root, 'describe', '--tags', '--always']).decode('utf-8').strip()
        
        # fetch current commit hash
        git_commit = subprocess.check_output(['git', '-C', repo_root, 'rev-parse', '--short', 'HEAD']).decode('utf-8').strip()
        
        # combine tag and commit 
        version = f"{git_tag}-{git_commit}"
        
        return version
    
    except subprocess.CalledProcessError:

        return "unknown"

# create the ROGUE1 PG header entry
def create_pg_header(header, command_line, description):
    
    version = get_version()
    
    # find the last PG entry to set as the previous program and determine the next ID
    last_pg = None
    next_id = 1
    if 'PG' in header:
        last_pg = header['PG'][-1]
        last_id_num = int(last_pg['ID'].split('.')[-1])
        next_id = last_id_num + 1

    new_header_entry = {
        'ID': f'ROGUE1.{next_id}',
        'PN': 'ROGUE1',
        'VN': version,
        'CL': command_line,
        'DS': description
    }
    
    if last_pg:
        new_header_entry['PP'] = last_pg['ID']
    
    return new_header_entry

def update_pg_header(header, command_line, description):
    new_header_entry = create_pg_header(header, command_line, description)
    
    if 'PG' in header:
        header['PG'].append(new_header_entry)
    else:
        header['PG'] = [new_header_entry]
    
    return header