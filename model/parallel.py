import sys
import subprocess

n_processes = 2
# Using https://stackoverflow.com/questions/19156467/
procs = []
print(f'Running {n_processes} processes of BeRCM in parallel...')
for i in range(n_processes):
    if sys.platform.startswith('win32'):
        proc = subprocess.Popen([sys.executable, 'run.py'], creationflags=subprocess.CREATE_NEW_CONSOLE)
    else:
        proc = subprocess.Popen([sys.executable, 'run.py'])
    procs.append(proc)

for proc in procs:
    proc.wait()

