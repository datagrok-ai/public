import gLPKConstructor from 'glpk.js';
import {GLPK, LP, Options, Result} from 'glpk.js';
import type {CobraModelData, CobraMetabolite, CobraReaction, CobraGene} from '../../escher_src/src/ts/types';

class GLPKWrapper {
  private static glpkInstance: GLPK;

  static async getInstance(): Promise<GLPK> {
    if (!GLPKWrapper.glpkInstance)
      GLPKWrapper.glpkInstance = await gLPKConstructor();

    return GLPKWrapper.glpkInstance;
  }
}

function prepareProblem(modelData: CobraModelData, glpk: GLPK): LP {
  const vars = modelData.reactions.map((reaction) => ({
    name: reaction.id,
    coef: reaction.objective_coefficient || 0,
    bnds: {
      type: reaction.lower_bound === reaction.upper_bound ? glpk.GLP_FX : glpk.GLP_DB,
      ub: reaction.upper_bound,
      lb: reaction.lower_bound
    }
  }));

  const subjectTo = modelData.metabolites.map((metabolite) => {
    const constraint = {
      name: metabolite.id,
      bnds: {type: glpk.GLP_FX, ub: 0, lb: 0}, // mass balance = 0
      vars: [] as Array<{name: string; coef: number}>
    };

    // Find all reactions that involve this metabolite
    modelData.reactions.forEach((reaction) => {
      if (reaction.metabolites[metabolite.id]) {
        constraint.vars.push({
          name: reaction.id,
          coef: reaction.metabolites[metabolite.id]
        });
      }
    });

    return constraint;
  });

  const lp: LP = {
    name: modelData.id || 'CobraModel',
    objective: {
      direction: glpk.GLP_MAX, // always maximize, if no positive coefficients, the objective will be zero
      name: 'obj',
      vars: vars.filter((v) => v.coef !== 0).map((v) => ({
        name: v.name,
        coef: v.coef
      }))
    },
    subjectTo,
    bounds: vars.map((v) => ({
      name: v.name,
      type: v.bnds.type,
      ub: v.bnds.ub,
      lb: v.bnds.lb
    }))
  };
  return lp;
}

export async function solveUsingGLPKJvail(modelData: CobraModelData): Promise<{fluxes: Float32Array[], reactionNames: string[]}> {
  const glpk = await GLPKWrapper.getInstance();

  // Build variables (reactions) with bounds
  const lp = prepareProblem(modelData, glpk);

  const options: Options = {
    msglev: glpk.GLP_MSG_OFF, // Disable output
    presol: true
  };
  // error on type side, needs to be async
  const sol = (await glpk.solve(lp, options)).result.vars;
  const reactionNames = modelData.reactions.map((r) => r.id);
  return {fluxes: [new Float32Array(reactionNames.map((rName) => sol[rName]))], reactionNames};
}

export async function extremeSolveUsingGLPKJvail(modelData: CobraModelData, start?: number, end?: number): Promise<{fluxes: Float32Array[], reactionNames: string[]}> {
  const glpk = await GLPKWrapper.getInstance();
  start ??= 0;
  end ??= modelData.reactions.length * 2;

  // Build variables (reactions) with bounds
  const lp = prepareProblem(modelData, glpk);

  const options: Options = {
    msglev: glpk.GLP_MSG_OFF, // Disable output
    presol: true,
  };

  const fluxes: Float32Array[] = [];
  let reactionNames: string[] = [];
  // set objective to each reaction in the range
  for (let i = start; i < end; i++) {
    const reactionIndex = Math.floor(i / 2);
    const isForward = (i % 2 === 0);
    const reaction = modelData.reactions[reactionIndex];
    // set objective
    lp.objective.vars = [{name: reaction.id, coef: isForward ? 1 : -1}];
    // solve, problem with types, need to be async
    const result = await glpk.solve(lp, options);
    // results.push(result.result.vars);
    const vs = result.result.vars;
    if (reactionNames.length == 0)
      reactionNames = modelData.reactions.map((r) => r.id);
    fluxes.push(new Float32Array(reactionNames.map((rName) => vs[rName])));
    // if (result?.result?.status !== glpk.GLP_OPT)
    //   throw new Error(`Optimization failed for reaction ${reaction.id} in ${isForward ? 'forward' : 'reverse'} direction.`);
  }
  return {fluxes, reactionNames};
}
