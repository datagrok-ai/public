#name: GlomEnvTestInplace
#environment: channels: [conda-forge], dependencies: [python=3.7, requests, dataclasses, glom, {pip: [requests]}]
#language: python
#output: string result

import json, os, base64, uuid, re, requests
import sys, getopt, datetime
from dataclasses import dataclass
from glom import glom
import pandas as pd

target = {'a': {'b': {'c': 'd'}}}
result = glom(target, 'a.b.c')  # returns 'd'