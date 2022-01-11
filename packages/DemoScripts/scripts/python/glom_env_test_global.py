#name: GlomEnvTest
#environment: DemoScripts:glom_global
#language: python
#top-menu: Test | Environment Test
#output: string result

import requests
from glom import glom
import pandas as pd

target = {'a': {'b': {'c': 'd'}}}
result = glom(target, 'a.b.c')  # returns 'd'