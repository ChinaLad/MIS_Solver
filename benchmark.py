import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Config
EXECUTABLE = "mis_v6"
INPUT = "small"
GRAPH_TYPE = ["sparse", "normal", "dense"]

time = 0
for t in GRAPH_TYPE:
    input_files = sorted(Path(f"data/{INPUT}/{t}").glob("*.txt"))
    for input_file in input_files:
        cmd = [f"build/{EXECUTABLE}", str(input_file)]
        print(f"Running: {' '.join(cmd)}")
        
        try:
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, timeout=120)
            output = output.decode("utf-8")

            # Parse output
            size_line = next(line for line in output.splitlines() if "MIS size:" in line)
            time_line = next(line for line in output.splitlines() if "Execution time:" in line)

            size = int(size_line.split(":")[1])
            time += float(time_line.split(":")[1].split()[0])

        except subprocess.TimeoutExpired:
            print(f"Timeout on {input_file}")
    
print(f"Results for ({EXECUTABLE}): {round(time/(len(GRAPH_TYPE)*10), 5)} seconds average")