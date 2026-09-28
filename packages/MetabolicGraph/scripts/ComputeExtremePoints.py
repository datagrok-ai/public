#name: ComputeExtremePoints
#description: Computes cobra's OptGP warmup points (minimum and maximum of every reaction flux) for the WebAssembly sampler
#language: python
#environment: channels: [Conda-forge], dependencies: [python=3.12, {pip: [cobra]}]
#input: string cobraModel
#output: string result
from cobra.io.dict import model_from_dict
from cobra.sampling import OptGPSampler
import json

model = model_from_dict(json.loads(cobraModel))

# cobra generates the warmup points while setting the sampler up. They live in the space of forward
# and reverse variables, which the WebAssembly sampler needs as is: each point lists the forward
# variables of all reactions, then the reverse ones.
sampler = OptGPSampler(model, processes=1, seed=42)
result = json.dumps({
    'reactionNames': [r.id for r in model.reactions],
    'points': [list(w[sampler.fwd_idx]) + list(w[sampler.rev_idx]) for w in sampler.warmup],
})
