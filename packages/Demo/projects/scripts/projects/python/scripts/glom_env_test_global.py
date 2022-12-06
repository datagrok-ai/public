#name: GlomEnvTest
#environment: Demo:Scripts:glom_global
#language: python
#output: string result

import requests
from glom import glom
import pandas as pd

target = {'a': {'b': {'c': 'd'}}}
result = glom(target, 'a.b.c')  # returns 'd'