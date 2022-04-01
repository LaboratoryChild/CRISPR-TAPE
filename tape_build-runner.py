import subprocess

if __name__ == '__main__':
    build_command = 'pyinstaller --onefile -c --icon=icon.ico -p . --python-option v -n CRISPR-TAPE CRISPR_TAPE/interface.py'
    subprocess.run(build_command, shell=True, check=True)