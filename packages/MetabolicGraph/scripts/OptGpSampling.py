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

# One chain, like the WebAssembly sampler, which returns these exact samples for the same seed.
# cobra's default runs one short chain per CPU core, which biases small runs towards the warmup points.
optgp = OptGPSampler(model, thinning, processes=1, seed=42)
res = optgp.sample(nSamples)
