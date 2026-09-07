#name: ComputeExtremePoints
#description: Computes extreme points (FBA for each reaction, maximize and minimize) to define the feasibility space
#language: python
#environment: channels: [Conda-forge], dependencies: [python=3.12, {pip: [cobra]}]
#input: string cobraModel
#output: string result
from cobra.io.dict import model_from_dict
from optlang.interface import OPTIMAL
import json

jsonMap = json.loads(cobraModel)
model = model_from_dict(jsonMap)

reaction_ids = [r.id for r in model.reactions]
n_reactions = len(reaction_ids)

# Build forward/reverse variable index mapping (same as cobra's HRSampler)
var_idx = {v: idx for idx, v in enumerate(model.variables)}
fwd_idx = [var_idx[r.forward_variable] for r in model.reactions]
rev_idx = [var_idx[r.reverse_variable] for r in model.reactions]

solutions = []
for sense in ('min', 'max'):
    model.objective_direction = sense
    for i, rxn in enumerate(model.reactions):
        fwd_var = model.variables[fwd_idx[i]]
        rev_var = model.variables[rev_idx[i]]

        # Skip fixed reactions
        if rxn.upper_bound - rxn.lower_bound < 1e-6:
            continue

        model.objective.set_linear_coefficients({fwd_var: 1, rev_var: -1})
        model.slim_optimize()

        if model.solver.status != OPTIMAL:
            model.objective.set_linear_coefficients({fwd_var: 0, rev_var: 0})
            continue

        # Extract net fluxes per reaction from solver primals (forward - reverse)
        primals = model.solver.primal_values
        fluxes = [
            primals[model.variables[fwd_idx[j]].name] - primals[model.variables[rev_idx[j]].name]
            for j in range(n_reactions)
        ]
        solutions.append(fluxes)

        # Reset objective coefficients
        model.objective.set_linear_coefficients({fwd_var: 0, rev_var: 0})

result = json.dumps({
    'reactionNames': reaction_ids,
    'solutions': solutions,
})
