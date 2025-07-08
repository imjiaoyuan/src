import os
import subprocess
from datetime import datetime

CONFIG = {
    'include_extensions': ['.md', '.ipynb', '.r', '.py', '.sh'],
    'exclude_files': ['update.py', 'README.md'],
    'exclude_dirs': ['.git', '__pycache__', 'src']
}

def normalize_path(path):
    return path.replace(os.path.sep, '/')

def generate_directory_tree(root_dir, base_dir, level=0):
    content = []
    indent = '  ' * level
    
    try:
        items = sorted(os.listdir(root_dir))
    except PermissionError:
        return content
    
    # Process directories first (to maintain proper tree structure)
    dirs = [d for d in items 
           if os.path.isdir(os.path.join(root_dir, d)) 
           and d not in CONFIG['exclude_dirs']]
    
    for dir_name in dirs:
        dir_path = os.path.join(root_dir, dir_name)
        sub_content = generate_directory_tree(dir_path, base_dir, level + 1)
        
        if sub_content:
            content.append(f"{indent}- **{dir_name}/**\n")
            content.extend(sub_content)
    
    # Then process files
    files = [f for f in items 
            if os.path.isfile(os.path.join(root_dir, f))
            and any(f.endswith(ext) for ext in CONFIG['include_extensions'])
            and f not in CONFIG['exclude_files']]
    
    for file in files:
        rel_path = os.path.relpath(os.path.join(root_dir, file), base_dir)
        content.append(f"{indent}- [{file}]({normalize_path(rel_path)})\n")
    
    return content

def git_commit():
    try:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        subprocess.run(["git", "add", "."], check=True)
        subprocess.run(["git", "commit", "-m", f"Update directory structure: {current_time}"], check=True)
        print("Git commit successful")
    except subprocess.CalledProcessError as e:
        print(f"Git error: {e}")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    
    content = ["# Repository Directory Structure\n\n"]
    content.extend(generate_directory_tree(parent_dir, parent_dir))
    
    readme_path = os.path.join(parent_dir, 'README.md')
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.writelines(content)
    
    print(f"Generated {readme_path}")
    git_commit()

if __name__ == "__main__":
    main()