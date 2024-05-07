import logging
import subprocess

# logging-enabled shell command function 
def run_command(cmd, message):
    logging.info(message)
    subprocess.run(cmd, shell=True, check=True)