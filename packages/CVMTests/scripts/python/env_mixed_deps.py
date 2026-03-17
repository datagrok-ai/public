#name: PythonEnvMixedDeps
#description: Environment with both conda extras and pip packages
#language: python
#environment: env_mixed_deps
#input: string data
#output: string result

import scipy
import numpy as np
import orjson
import msgpack
from scipy.stats import norm

# Verify conda dep (scipy) works with numpy
mean, std = 0, 1
cdf_val = float(norm.cdf(1.96, mean, std))

# Verify pip dep (orjson) works
json_bytes = orjson.dumps({'key': data, 'cdf': round(cdf_val, 4)})
parsed = orjson.loads(json_bytes)

# Verify pip dep (msgpack) works
packed = msgpack.packb({'test': data})
unpacked = msgpack.unpackb(packed)

result = parsed['key'] + ':' + str(parsed['cdf']) + ':' + unpacked['test']
