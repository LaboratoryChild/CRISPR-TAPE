import os
import shutil
import subprocess
import sys

def create_desktop_file(executable_path, icon_path):
    """
    Creates a .desktop file for the application with the specified icon.

    executable_path: The full path to the executable.
    icon_path: The full path to the icon file.
    """
    desktop_entry = f"""[Desktop Entry]
                    Version=2.0
                    Type=Application
                    Name=CRISPR-TAPE
                    Exec={executable_path}
                    Icon={icon_path}
                    Terminal=false
                    Categories=Utility;
                    """
    desktop_file_path = os.path.expanduser('~/.local/share/applications/CRISPR-TAPE.desktop')
    with open(desktop_file_path, 'w') as desktop_file:
        desktop_file.write(desktop_entry)
    print(f'Created .desktop file at {desktop_file_path}')

def build_executable(script_dir, icon_path):
    """
    Builds the executable using PyInstaller.

    script_dir: The directory containing the CRISPR-TAPE python files.
    icon_path: The path to the icon file.
    """
    # Prepare the PyInstaller command
    pyinstaller_command = f'pyinstaller -w --onefile --paths venv/lib/site-packages --icon={icon_path} -p . '
    pyinstaller_command += f'-n CRISPR-TAPE {os.path.join(script_dir, "__main__.py")} {os.path.join(script_dir, "general_function.py")} '
    pyinstaller_command += f'{os.path.join(script_dir, "shared_functions.py")} {os.path.join(script_dir, "specific_function.py")} {os.path.join(script_dir, "__init__.py")}'
    # Use PyInstaller to build the executable
    subprocess.run(pyinstaller_command, shell=True, check=True)
    # move the executable to the current directory
    if os.path.exists("CRISPR-TAPE"):
        os.remove("CRISPR-TAPE")
    if os.path.exists("dist/CRISPR-TAPE"):
        shutil.move("dist/CRISPR-TAPE", ".")

def main():
    """
    Main function to clone a repo and build an executable from it.
    """
    # Define the directory where the repo should be cloned
    script_dir = 'CRISPR_TAPE'  # Replace with the path to your main script
    icon_path = 'icon.ico'  # Replace with your actual icon path
    assert os.path.exists(icon_path) and os.path.exists(script_dir)
    # use the absolute paths
    script_dir = os.path.abspath(script_dir)
    icon_path = os.path.abspath(icon_path)
    print('Building executable with icon...')
    build_executable(script_dir, icon_path)
    # make a .desktop file if this is a linux system
    if sys.platform == "linux":
        create_desktop_file(os.path.realpath("CRISPR-TAPE"), os.path.realpath(icon_path))
    print('Build completed.')

if __name__=="__main__":
    main()
