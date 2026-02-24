
import {
  glp_create_prob, glp_set_prob_name, glp_set_obj_dir, glp_add_rows,
  glp_add_cols, glp_set_row_name, glp_set_row_bnds, glp_set_col_name,
  glp_set_col_bnds, glp_set_obj_coef, glp_load_matrix, glp_simplex,
  glp_get_obj_val, glp_get_num_cols, glp_get_col_name, glp_get_col_prim,
  SMCP, GLP_MAX, GLP_FX, GLP_DB, GLP_ON
} from 'glpk.js'
import {Model, Solution} from './cobra-model';
export {Model, Solution, modelFromJsonData, modelFromWorkerData} from './cobra-model';

export function buildGlpkProblem (model: Model) {
  const nRows = model.metabolites.length
  const nCols = model.reactions.length
  const ia: number[] = []
  const ja: number[] = []
  const ar: number[] = []
  const metLookup: Record<string, number> = {}

  var lp = glp_create_prob()
  glp_set_prob_name(lp, 'knockout FBA')
  glp_set_obj_dir(lp, GLP_MAX)
  glp_add_rows(lp, nRows)
  glp_add_cols(lp, nCols)

  model.metabolites.forEach(function (metabolite, i) {
    var rowInd = i + 1
    glp_set_row_name(lp, rowInd, metabolite.id)
    glp_set_row_bnds(lp, rowInd, GLP_FX, 0.0, 0.0)
    metLookup[metabolite.id] = rowInd
  })

  var matInd = 1
  model.reactions.forEach(function (reaction, i) {
    var colInd = i + 1

    glp_set_col_name(lp, colInd, reaction.id)
    if (reaction.lower_bound === reaction.upper_bound)
      glp_set_col_bnds(lp, colInd, GLP_FX, reaction.lower_bound, reaction.upper_bound)
    else
      glp_set_col_bnds(lp, colInd, GLP_DB, reaction.lower_bound, reaction.upper_bound)

    if ('objective_coefficient' in reaction && !!reaction.objective_coefficient)
      glp_set_obj_coef(lp, colInd, reaction.objective_coefficient)

    for (var metId in reaction.metabolites) {
      ia[matInd] = metLookup[metId]
      ja[matInd] = colInd
      ar[matInd] = reaction.metabolites[metId]
      matInd++
    }
  })
  glp_load_matrix(lp, ia.length - 1, ia, ja, ar)

  return lp
}

export function optimizeModel (model: Model) {
  const problem = buildGlpkProblem(model)
  var smcp = new SMCP({ presolve: GLP_ON })
  const returnCode = glp_simplex(problem, smcp)
  var f = null
  var x = null
  if (returnCode === 0) {
    f = glp_get_obj_val(problem)
    x = {} as Record<string, number>
    for (var i = 1; i <= glp_get_num_cols(problem); i++)
      x[glp_get_col_name(problem, i)] = glp_get_col_prim(problem, i)
  } else {
    throw new Error(`INVALID SOLUTION. GLPK optimization failed with code ${returnCode}`)
  }

  return new Solution(f, x)
}

export function sampleExtremePoints (model: Model) {
  const problem = buildGlpkProblem(model);
  var smcp = new SMCP({ presolve: GLP_ON });
  const solutions: Solution[] = [];
  const errors: number[] = [];
  for (let i = 0; i < model.reactions.length * 2; i++) {
    const reactionIndex = Math.floor(i / 2);
    const direction = (i % 2 === 0) ? 1 : -1;
    for (let j = 1; j <= model.reactions.length; j++)
      glp_set_obj_coef(problem, j, 0);
    glp_set_obj_coef(problem, reactionIndex + 1, direction);
    const returnCode = glp_simplex(problem, smcp);
    let f = null;
    let x = null;
    if (returnCode === 0) {
      f = glp_get_obj_val(problem)
      x = {} as Record<string, number>
      for (let k = 1; k <= glp_get_num_cols(problem); k++)
        x[glp_get_col_name(problem, k)] = glp_get_col_prim(problem, k)
      solutions.push(new Solution(f, x));
    } else {
      errors.push(returnCode);
    }
  }
  return {solutions, errors};
}
