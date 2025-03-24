#!/usr/bin/env python
# This script merges multiple YAML CloudFormation template files into one,
# optionally saving the result to a user-defined file or auto-generating an output name.

import argparse
import yaml
from deepmerge import always_merger
import os

parser = argparse.ArgumentParser()
# Accepts one or more input YAML template files
parser.add_argument(
    "-i",
    "--input",
    help="Path to the template file",
    default=[], action='append'
)
# Optional output file path
parser.add_argument(
    "-o",
    "--output",
    help="Path for the result file",
    default=''
)
args = parser.parse_args()

result = {}
output_names=[]
# Iterate through each provided input file
for i in args.input:
    with open(i, 'r') as file:
        data = yaml.safe_load(file)
    result = always_merger.merge(result, data)
    output_names += os.path.basename(i).replace('.template.yml', '')

if args.output:
    output = args.output
else:
    output = "-".join(output_names) + '.yml'

# Write the merged YAML to the output file
with open(args.output, 'w') as file:
    yaml.dump(result, file, explicit_start=True)
