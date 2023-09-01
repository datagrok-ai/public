#!/usr/bin/env python

import argparse
import yaml
from deepmerge import always_merger
import os

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    help="Path to the template file",
    default=[], action='append'
)
parser.add_argument(
    "-o",
    "--output",
    help="Path for the result file",
    default=''
)
args = parser.parse_args()

result = {}
output_names=[]
for i in args.input:
    with open(i, 'r') as file:
        data = yaml.safe_load(file)
    result = always_merger.merge(result, data)
    output_names += os.path.basename(i).replace('.template.yml', '')

if args.output:
    output = args.output
else:
    output = "-".join(output_names) + '.yml'

with open(args.output, 'w') as file:
    yaml.dump(result, file, explicit_start=True)
