#!/usr/bin/env python3
import os, sys, boto3, logging, argparse
from concurrent.futures import ThreadPoolExecutor

LOCAL_DIRECTORY_PATH = '/home/jy/work/notebook'
OTHERS_NOTES_PATH = '/home/jy/work/notebook_others'
USERNAME = 'testuser'

try:
    with open(os.path.join(os.path.dirname(__file__), '.env'), 'r') as f:
        config = {k.strip(): v.strip().strip("'\"") for k, v in (line.split('=', 1) for line in f if '=' in line)}
except FileNotFoundError:
    print("Error: '.env' file not found.")
    sys.exit(1)

R2_BUCKET_NAME = config.get('R2_BUCKET_NAME')
R2_ACCOUNT_ID = config.get('R2_ACCOUNT_ID')
R2_ACCESS_KEY_ID = config.get('R2_ACCESS_KEY_ID')
R2_SECRET_ACCESS_KEY = config.get('R2_SECRET_ACCESS_KEY')
MAX_WORKERS = int(config.get('MAX_WORKERS', 10))

if not all([R2_BUCKET_NAME, R2_ACCOUNT_ID, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY]):
    print("Error: Required R2 variables are missing from '.env' file.")
    sys.exit(1)

logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(message)s', 
                    filename=os.path.join(os.path.dirname(__file__), 'log'), filemode='a')

def get_local_map(path):
    if not os.path.exists(path): os.makedirs(path)
    return {os.path.relpath(p, path).replace(os.path.sep, '/'): {'path': p, 'size': os.path.getsize(p)}
            for r, _, f in os.walk(path) for file in f if (p := os.path.join(r, file))}

def execute_sync(action):
    try:
        s3 = boto3.client('s3', endpoint_url=f"https://{R2_ACCOUNT_ID}.r2.cloudflarestorage.com",
                          aws_access_key_id=R2_ACCESS_KEY_ID, aws_secret_access_key=R2_SECRET_ACCESS_KEY, region_name="auto")
        s3.head_bucket(Bucket=R2_BUCKET_NAME)
    except Exception as e:
        logging.error(f"Initialization failed: {e}"); return False

    is_successful = True
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []
        if action == 'push':
            local_map = get_local_map(LOCAL_DIRECTORY_PATH)
            prefix = f"{USERNAME}/"
            paginator = s3.get_paginator('list_objects_v2')
            r2_map = {obj['Key']: {'size': obj['Size']} for page in paginator.paginate(Bucket=R2_BUCKET_NAME, Prefix=prefix) for obj in page.get('Contents', [])}
            
            push_local_map = {f"{prefix}{k}": v for k, v in local_map.items()}
            keys_to_delete = list(r2_map.keys() - push_local_map.keys())
            if keys_to_delete:
                for i in range(0, len(keys_to_delete), 1000):
                    executor.submit(s3.delete_objects, Bucket=R2_BUCKET_NAME, Delete={'Objects': [{'Key': k} for k in keys_to_delete[i:i+1000]]})
            
            keys_to_upload = [k for k, v in push_local_map.items() if k not in r2_map or v['size'] != r2_map.get(k,{}).get('size')]
            for key in keys_to_upload: futures.append(executor.submit(s3.upload_file, push_local_map[key]['path'], R2_BUCKET_NAME, key))
        else:
            paginator = s3.get_paginator('list_objects_v2')
            full_r2_map = {obj['Key']: {'size': obj['Size']} for page in paginator.paginate(Bucket=R2_BUCKET_NAME) for obj in page.get('Contents', [])}

            my_notes_r2 = {k.split('/', 1)[1]: v for k, v in full_r2_map.items() if k.startswith(f"{USERNAME}/")}
            my_notes_local = get_local_map(LOCAL_DIRECTORY_PATH)
            for key in (my_notes_local.keys() - my_notes_r2.keys()):
                try: os.remove(my_notes_local[key]['path'])
                except OSError as e: logging.error(f"Local delete failed for {key}: {e}"); is_successful = False
            for key, data in my_notes_r2.items():
                if key not in my_notes_local or data['size'] != my_notes_local.get(key, {}).get('size'):
                    local_path = os.path.join(LOCAL_DIRECTORY_PATH, key.replace('/', os.path.sep))
                    os.makedirs(os.path.dirname(local_path), exist_ok=True)
                    futures.append(executor.submit(s3.download_file, R2_BUCKET_NAME, f"{USERNAME}/{key}", local_path))
            
            others_notes_r2 = {k: v for k, v in full_r2_map.items() if not k.startswith(f"{USERNAME}/") and '/' in k}
            others_notes_local = get_local_map(OTHERS_NOTES_PATH)
            for key in (others_notes_local.keys() - others_notes_r2.keys()):
                try: os.remove(others_notes_local[key]['path'])
                except OSError as e: logging.error(f"Local delete failed for {key}: {e}"); is_successful = False
            for key, data in others_notes_r2.items():
                if key not in others_notes_local or data['size'] != others_notes_local.get(key, {}).get('size'):
                    local_path = os.path.join(OTHERS_NOTES_PATH, key.replace('/', os.path.sep))
                    os.makedirs(os.path.dirname(local_path), exist_ok=True)
                    futures.append(executor.submit(s3.download_file, R2_BUCKET_NAME, key, local_path))

        for future in futures:
            try: future.result()
            except Exception as e: logging.error(f"A sync task failed: {e}"); is_successful = False
    return is_successful

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['push', 'pull'])
    args = parser.parse_args()
    
    if execute_sync(args.action):
        print("Sync operation successful.")
    else:
        print("Sync operation failed. See log for details.")
        sys.exit(1)