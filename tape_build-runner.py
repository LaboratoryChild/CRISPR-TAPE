import subprocess

if __name__ == '__main__':
    build_command = 'pyinstaller --onefile -w -n CRISPR-TAPE CRISPR_TAPE/interface.py'
    subprocess.run(build_command, shell=True, check=True)