#name: OptGpSampling
#description: Samples the metabolic map and returns raw flux samples
#language: python
#environment: channels: [Conda-forge], dependencies: [python=3.12, {pip: [cobra]}]
#input: string cobraModel
#input: int nSamples = 1000 {nullable: true}
#input: int thinning = 1 {nullable: true}
#output: dataframe res
from cobra.sampling import OptGPSampler
from cobra.io.dict import model_from_dict
import json

jsonMap = json.loads(cobraModel)
model = model_from_dict(jsonMap)

optgp = OptGPSampler(model, thinning, seed=42, n_samples=nSamples)
res = optgp.sample(nSamples)