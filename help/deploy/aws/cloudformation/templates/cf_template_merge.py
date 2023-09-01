#!/usr/bin/env python

import argparse
import yaml
from deepmerge import always_merger

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
    default='cloudformation.yaml'
)
args = parser.parse_args()

result = {}
for i in args.input:
    with open(i, 'r') as file:
        data = yaml.safe_load(file)
    result = always_merger.merge(result, data)

with open(args.output, 'w') as file:
    yaml.dump(result, file, explicit_start=True)
