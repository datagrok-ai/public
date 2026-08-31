/* eslint-disable camelcase */
/* eslint-disable new-cap */
import glpkWasm from 'glpk-wasm';
import type {CobraModelData, CobraMetabolite, CobraReaction, CobraGene} from '../../../escher_src/src/ts/types';

class GLPKSingleton {
  private static instance: GLPKSingleton;
  glpk: ReturnType<Awaited<typeof glpkWasm>>;
  private constructor() { }
  static getInstance() {
    if (!GLPKSingleton.instance) {
      GLPKSingleton.instance = new GLPKSingleton();
      // @ts-ignore
      GLPKSingleton.instance.glpk = glpkWasm({locateFile: () => 'glpk.all.wasm'});
    }
    return GLPKSingleton.instance;
  }
}

export async function runGLPKFBA(network: CobraModelData, findSolutionSpaceEdges = false, start?: number, end?: number) {
  start ??= 0;
  end ??= network.reactions.length * 2;
  const glpk = await GLPKSingleton.getInstance().glpk;
  const GLP_SF_AUTO = 128;
  const {_glp_create_prob, _glp_set_prob_name, _glp_set_obj_dir, _glp_add_rows,
    _glp_add_cols, _glp_set_row_name, _glp_set_row_bnds, _glp_set_col_name,
    _glp_set_col_bnds, _glp_set_obj_coef, _glp_get_obj_coef, _glp_load_matrix, _glp_simplex,
    _glp_get_obj_val, _glp_get_num_cols, _glp_get_col_name, _glp_get_col_prim,
    _glp_init_smcp, _glp_scale_prob, _free, _malloc, getValue, setValue, stringToUTF8, UTF8ToString, lengthBytesUTF8, _glp_copy_prob, _glp_erase_prob, _glp_delete_prob} = glpk;
  glpk._glp_term_out(0); // turn off glpk output
  class Model {
    id: string = '';
    name: string = '';
    description: string = '';
    notes: any;
    metabolites: CobraMetabolite[] = [];
    reactions: CobraReaction[] = [];
    genes: CobraGene[] = [];
    // constructor (data) {
    //     this.reactions = data.reactions.map(x => ({...x}))
    //     this.metabolites = data.metabolites.map(x => ({...x}))
    //     this.genes = data.genes.map(x => ({...x}))
    //     this.id = data.id
    //     this.notes = data.notes // TODO is this an object? if so clone
    //     this.description = data.description
    // }

    buildGlpkProblem() {
      /** Build a GLPK LP for the model */
      //console.time('build LP')
      const nRows = this.metabolites.length;
      const nCols = this.reactions.length;
      let ia: any = [];
      let ja: any = [];
      let ar: any = [];
      const metLookup = {};

      // initialize LP objective
      const lp = _glp_create_prob();
      _glp_set_prob_name(lp, getStringPointer('knockout FBA'));
      // maximize
      _glp_set_obj_dir(lp, GLP_MAX);
      // set up rows and columns
      _glp_add_rows(lp, nRows);
      _glp_add_cols(lp, nCols);

      // metabolites
      this.metabolites.forEach(function(metabolite, i) {
        const rowInd = i + 1;
        _glp_set_row_name(lp, rowInd, getStringPointer(metabolite.id));
        _glp_set_row_bnds(lp, rowInd, GLP_FX, 0.0, 0.0);
        // remember the indices of the metabolites
        metLookup[metabolite.id] = rowInd;
      });

      // reactions
      let matInd = 1;
      this.reactions.forEach(function(reaction, i) {
        const colInd = i + 1;

        _glp_set_col_name(lp, colInd, getStringPointer(reaction.id));
        if (reaction.lower_bound === reaction.upper_bound)
          _glp_set_col_bnds(lp, colInd, GLP_FX, reaction.lower_bound, reaction.upper_bound);
        else
          _glp_set_col_bnds(lp, colInd, GLP_DB, reaction.lower_bound, reaction.upper_bound);


        // object_coefficient is optional for reaction in COBRA JSON
        if ('objective_coefficient' in reaction && !!reaction.objective_coefficient) {
          //console.log(reaction.id, reaction.objective_coefficient);
          _glp_set_obj_coef(lp, colInd, reaction.objective_coefficient);
        }

        // S matrix values
        for (const metId in reaction.metabolites) {
          ia[matInd] = metLookup[metId];
          ja[matInd] = colInd;
          ar[matInd] = reaction.metabolites[metId];
          matInd++;
        }
      });
      //console.log(matInd)
      // Load the S matrix
      // console.log(lp);
      ia = Int32Array.from(ia);
      ja = Int32Array.from(ja);
      ar = Float64Array.from(ar);
      const iaPointer = mallocWithStorage(ia.length * ia.BYTES_PER_ELEMENT);
      const jaPointer = mallocWithStorage(ja.length * ja.BYTES_PER_ELEMENT);
      const arPointer = mallocWithStorage(ar.length * ar.BYTES_PER_ELEMENT);
      for (let i = 0; i < ia.length; i++) {
        setValue(iaPointer + i * ia.BYTES_PER_ELEMENT, ia[i], 'i32');
        setValue(jaPointer + i * ja.BYTES_PER_ELEMENT, ja[i], 'i32');
        setValue(arPointer + i * ar.BYTES_PER_ELEMENT, ar[i], 'double');
      }

      // setValue(iaPointer + 4, 12, 'i32');

      // console.log(getValue(iaPointer + 4));
      //console.log(ia.length - 1, ia, ja, ar);
      _glp_load_matrix(lp, ia.length - 1, iaPointer, jaPointer, arPointer);
      //console.log(_glp_load_matrix.prototype);
      // console.time('glp_scale_prob')

      // console.timeEnd('glp_scale_prob')
      //console.timeEnd('build LP')
      return lp;
    }

    // here reaction_number is 0 based index
    set_lone_objective(problem, reaction_number, direction) {
      this.reactions.forEach((r, i) => {
        _glp_set_obj_coef(problem, i + 1, 0);
      });
      _glp_set_obj_coef(problem, reaction_number + 1, direction);
    }

    /** Pass the problem from build step */
    optimize(problem) {
      _glp_scale_prob(problem, GLP_SF_AUTO);
      //console.time('optimize LP')
      // const problem = this.buildGlpkProblem();
      // Use default simplex parameters by passing NULL pointer
      //console.log(smcp);
      //console.timeEnd('optimize LP')
      // console.log('coef ', _glp_get_obj_coef(problem, 1));
      //console.time('optimize LP')
      const returnCode = _glp_simplex(problem, 0 as unknown as any);


      //console.timeEnd('optimize LP')
      // _glp_set_obj_coef(problem, 1, -1);
      // console.log('coef ', _glp_get_obj_coef(problem, 1));
      let f = null;
      let x = null;
      if (returnCode === 0) {
        // get the objective
        f = _glp_get_obj_val(problem);
        console.log('obj val ', f);
        // get the primal
        x = {};
        //console.log( _glp_get_num_cols(problem));
        for (let i = 1; i <= _glp_get_num_cols(problem); i++) {
          //console.log(i, UTF8ToString(_glp_get_col_name(problem, i), 256))
          x[UTF8ToString(_glp_get_col_name(problem, i))] = _glp_get_col_prim(problem, i);
        }
      } else
        throw new Error(`INVALID SOLUTION. GLPK optimization failed with code ${returnCode}`);


      return new Solution(f, x);
    }

    slim_optimize(problem) {
      _glp_scale_prob(problem, GLP_SF_AUTO);
      const returnCode = _glp_simplex(problem, 0 as unknown as any);
      if (returnCode === 0) {
        const fluxes = new Float32Array(this.reactions.length);
        for (let i = 1; i <= _glp_get_num_cols(problem); i++)
          fluxes[i - 1] = _glp_get_col_prim(problem, i);

        return fluxes;
      } else
        throw new Error(`INVALID SOLUTION. GLPK optimization failed with code ${returnCode}`);
    }
  }

  class Solution {
    fluxes: { [key: string]: number; };
    objectiveValue: number;
    constructor(objectiveValue, fluxes) {
      this.objectiveValue = objectiveValue;
      this.fluxes = fluxes;
    }
  }

  const GLP_DB = 4;
  const GLP_FX = 5;
  const GLP_MAX = 2;
  const GLP_ON = 1;

  const pointerStorage = [];

  function mallocWithStorage(size: number): number {
    const ptr = _malloc(size);
    pointerStorage.push(ptr);
    return ptr as number;
  }

  function getStringPointer(str: string) {
    const ptr = mallocWithStorage(lengthBytesUTF8(str) + 1);
    stringToUTF8(str, ptr, lengthBytesUTF8(str) + 1);
    return ptr;
  }

  function modelFromJsonData(data: CobraModelData): Model | null {
    if (data === null) return null;

    const model = new Model();

    // when objective coefficients are not -1, 1, or 0, change them to 1/-1
    model.reactions = data.reactions.map((r) => {
      let coeff = 0;
      if (r.objective_coefficient && r.objective_coefficient !== 0) {
        if (r.objective_coefficient < 0) coeff = -1;
        else coeff = 1;
      }
      return ({...r, objective_coefficient: coeff});
    });

    model.metabolites = data.metabolites.map((x) => ({...x}));
    model.genes = data.genes.map((x) => ({...x}));
    model.id = data.id;
    model.notes = data.notes; // TODO is this an object? if so clone
    model.description = data.description;
    return model;
  }

  const lbs = new Float32Array(network.reactions.map((r) => r.lower_bound));
  const ubs = new Float32Array(network.reactions.map((r) => r.upper_bound));
  const model = modelFromJsonData(network);
  const problem = model.buildGlpkProblem();


  const solutions: Float32Array[] = [];
  let reactionNames: string[] | null = null;

  console.time('optimizations');
  if (findSolutionSpaceEdges) {
    for (let i = start; i < end; i += 1) {
      const newProblemPtr = problem;
      //_glp_copy_prob(newProblemPtr, problem, GLP_ON);
      // _glp_erase_prob(problem);
      //const model = modelFromJsonData(network);
      const j = Math.floor(i / 2);
      const direction = (i % 2 === 0) ? 1 : -1;
      model.set_lone_objective(newProblemPtr, j, direction);
      // set all objective coefficients to 0
      //model.reactions.forEach(r => r.objective_coefficient = 0);
      // set the j-th reaction objective coefficient to direction
      //model.reactions[j].objective_coefficient = direction;
      //const problem = model.buildGlpkProblem();
      const solution = model.slim_optimize(newProblemPtr);
      if (!reactionNames)
        reactionNames = network.reactions.map((r) => r.id);
      // console.log(solution);
      // break;
      solutions.push(solution);
      if (i % 50 === 0)

        console.log(`completed ${(i / network.reactions.length / 2 * 100).toFixed(4)}% optimizations`);

      // _glp_delete_prob(newProblemPtr);
    }
  } else {
    const solution = model.optimize(problem).fluxes;
    solutions.push(new Float32Array(Object.values(solution)));
    reactionNames = Object.keys(solution);
  }


  console.timeEnd('optimizations');
  pointerStorage.forEach((p) => _free(p));
  return {lbs, ubs, solutions, reactionNames};
}
