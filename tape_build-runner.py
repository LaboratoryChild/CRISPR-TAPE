import subprocess

if __name__ == '__main__':
    build_command = 'pyinstaller -w --paths venv/lib/site-packages --icon=icon.ico -p . -n CRISPR-TAPE CRISPR_TAPE/interface.py CRISPR_TAPE/shared_functions.py CRISPR_TAPE/general_function.py CRISPR_TAPE/specific_function.py'
    subprocess.run(build_command, shell=True, check=True)