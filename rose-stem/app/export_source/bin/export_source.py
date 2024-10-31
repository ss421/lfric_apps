#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

import os
import subprocess

def run_command(command):
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    """
    result = subprocess.run(
        command.split(),
        capture_output=True,
        text=True,
        timeout=600,
        shell=False,
        check=False,
    )
    print(result.stdout)
    if result.returncode:
        raise RuntimeError(
            f"The command '{command}' failed with error:\n\n{result.stderr}"
        )

# Read paths and sources from environment variables
source_root = os.environ["SOURCE_ROOT"]
apps_source = os.environ["SOURCE_LFRIC_APPS"]
apps_path = os.path.join(source_root, "apps")
core_source = os.environ["SOURCE_LFRIC_CORE"]
core_path = os.path.join(source_root, "core")
ctldata_path = os.path.join(source_root, "ctldata")

print(f"\n[INFO] LFRic Apps Source: {apps_source}")
print(f"[INFO] LFRic Core Source: {core_source}\n")

# Create directories for the different sources
mkdir_command = f"mkdir -p "
run_command(f"{mkdir_command}{apps_path}")
run_command(f"{mkdir_command}{core_path}")
run_command(f"{mkdir_command}{ctldata_path}")

# Export entire apps and core directories
run_command(f"fcm export --force {apps_source} {apps_path}")
run_command(f"fcm export --force {core_source} {core_path}")

# Read through the file defining ctldata sources
# For each source, export into ctldata directory
dests = []
ctldata_file = "ctldata_list.txt"
with open(ctldata_file) as f:
    for line in f:
        line = line.strip()
        line = line.strip("\n")
        if line.startswith("#") or len(line) == 0:
            continue
        source, dest = line.split()
        if dest in dests:
            raise RuntimeError(
                "Two sources in the ctldata_list.txt file share the same "
                f"destination: {dest}. Please deduplicate these."
            )
        dests.append(dest)
        print(f"[INFO] Source: {source}   Dest: {dest}")
        run_command(
            f"fcm export --force {source} {os.path.join(ctldata_path, dest)}"
            )

