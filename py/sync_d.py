import os
import shutil
import datetime
import filecmp

def sync_path(from_path, to_path):
    changes = []

    if os.path.isfile(from_path):
        to_dir = os.path.dirname(to_path)
        if not os.path.exists(to_dir):
            os.makedirs(to_dir)
            changes.append(f"CREATED DIR: {to_dir}")
        
        if not os.path.exists(to_path) or not filecmp.cmp(from_path, to_path, shallow=True):
            shutil.copy2(from_path, to_path)
            changes.append(f"UPDATED FILE: {to_path}")
        return changes

    if not os.path.isdir(from_path):
        return changes

    if not os.path.exists(to_path):
        os.makedirs(to_path)
        changes.append(f"CREATED DIR: {to_path}")

    for src_dir, _, files in os.walk(from_path):
        dest_dir = src_dir.replace(from_path, to_path, 1)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
            changes.append(f"CREATED DIR: {dest_dir}")

        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dest_file = os.path.join(dest_dir, file_)
            if not os.path.exists(dest_file) or not filecmp.cmp(src_file, dest_file, shallow=True):
                shutil.copy2(src_file, dest_file)
                changes.append(f"UPDATED FILE: {dest_file}")

    for dest_dir, dirs, files in os.walk(to_path, topdown=False):
        for file_ in files:
            dest_file = os.path.join(dest_dir, file_)
            src_file = dest_file.replace(to_path, from_path, 1)
            if not os.path.exists(src_file):
                os.remove(dest_file)
                changes.append(f"DELETED FILE: {dest_file}")
        
        for dir_ in dirs:
            dest_subdir = os.path.join(dest_dir, dir_)
            src_subdir = dest_subdir.replace(to_path, from_path, 1)
            if not os.path.exists(src_subdir):
                shutil.rmtree(dest_subdir)
                changes.append(f"DELETED DIR: {dest_subdir}")
                
    return changes

if __name__ == "__main__":
    backup_paths = [
        {
            "from": "/mnt/c/Users/JiaoYuan/Documents/NetSarang Computer",
            "to": "/mnt/d/BACKUP/JiaoYuan/Documents/NetSarang Computer"
        },
        {
            "from": "/mnt/c/Users/JiaoYuan/Documents/My Games",
            "to": "/mnt/d/BACKUP/JiaoYuan/Documents/My Games"
        },
        {
            "from": "/mnt/c/Users/JiaoYuan/Pictures",
            "to": "/mnt/d/BACKUP/JiaoYuan/Pictures"
        },
        {
            "from": "/mnt/c/Users/JiaoYuan/.ssh",
            "to": "/mnt/d/BACKUP/JiaoYuan/.ssh"
        },
        {
            "from": "/mnt/c/Users/JiaoYuan/.wslconfig",
            "to": "/mnt/d/BACKUP/JiaoYuan/.wslconfig"
        },
        {
            "from": "/mnt/c/User/JiaoYuan/Desktop/MISC/CV.pdf",
            "to": "/mnt/d/BACKUP/JiaoYuan/Desktop/MISC/CV.pdf"
        },
        {
            "from": "/mnt/c/Scoop/apps/rclone/current/rclone.conf",
            "to": "/mnt/d/BACKUP/Scoop/apps/rclone/rclone.conf"
        },
        {
            "from": "/mnt/c/Scoop/apps/vscode/current/data/user-data/User/settings.json",
            "to": "/mnt/d/BACKUP/Scoop/apps/vscode/current/data/user-data/User/settings.json"
        },
        {
            "from": "/home/jy/.bashrc",
            "to": "/mnt/d/BACKUP/jy/.bashrc"
        },
        {
            "from": "/home/jy/.gitconfig",
            "to": "/mnt/d/BACKUP/jy/.gitconfig"
        },
        {
            "from": "/home/jy/.git-credentials",
            "to": "/mnt/d/BACKUP/jy/.git-credentials"
        },
        {
            "from": "/home/jy/.condarc",
            "to": "/mnt/d/BACKUP/jy/.condarc"
        },
        {
            "from": "/home/jy/.Rprofile",
            "to": "/mnt/d/BACKUP/jy/.Rprofile"
        },
        {
            "from": "/home/jy/.npmrc",
            "to": "/mnt/d/BACKUP/jy/.npmrc"
        },
        {
            "from": "/home/jy/.config/rclone",
            "to": "/mnt/d/BACKUP/jy/.config/rclone"
        },
        {
            "from": "/home/jy/work/2fa",
            "to": "/mnt/d/BACKUP/jy/work/2fa"
        }
    ]

    print(f"Backup started at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    total_changes = 0

    for paths in backup_paths:
        from_path = paths["from"]
        to_path = paths["to"]

        if not os.path.exists(from_path):
            print(f"WARNING: Source path does not exist, skipping: {from_path}")
            continue
        
        try:
            changes_made = sync_path(from_path, to_path)
            if changes_made:
                total_changes += len(changes_made)
                print(f"\n[SYNC] '{from_path}' -> '{to_path}'")
                for change in changes_made:
                    print(f"  - {change}")
        except Exception as e:
            print(f"ERROR: A critical error occurred while syncing '{from_path}': {e}")


    print(f"\nBackup finished. Total changes: {total_changes}")