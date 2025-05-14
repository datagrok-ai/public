/* eslint-disable valid-jsdoc */
// The MOEA/D method for multi-objective optimization: https://ieeexplore.ieee.org/document/4358754

//@ts-ignore: no types
import * as jStat from 'jstat';
import {DEFAULT_SETTINGS, Func, InputOptions, MoeadOptions, MoeadOutput} from './defs';
import {clip, euclideanDistance, pickTwo} from './utils';

export class Moead {
  public objective: Func;
  private inputOpts: InputOptions;
  private outputDim: number;
  private methodOpts: MoeadOptions;

  constructor(
    objective: Func,
    inputOpts: InputOptions,
    outputDim: number,
    methodOpts: MoeadOptions = DEFAULT_SETTINGS) {
    this.objective = objective;
    this.inputOpts = inputOpts;
    this.outputDim = outputDim;
    this.methodOpts = methodOpts;
  }

  /** Scalarization: weighted sum */
  public weightedSum(x: Float32Array, w: Float32Array): number {
    let sum = 0;
    const out = this.objective(x);

    for (let i = 0; i < this.outputDim; ++i)
      sum += w[i] * out[i];

    return sum;
  }

  /** Generate uniform weight vectors */
  public generateWeights(): Float32Array[] {
    const nWeights = this.methodOpts.nWeights;
    const weights = new Array<Float32Array>(nWeights + 1);

    for (let i = 0; i <= nWeights; ++i) {
      const w = new Float32Array(this.outputDim);

      let rem = 1;

      for (let k = 0; k < this.outputDim - 1; ++k) {
        w[k] = i / nWeights;
        rem -= w[k];
      }

      w[this.outputDim - 1] = rem;

      weights[i] = w;
    }

    return weights;
  }

  /**  Initialize a random solution */
  public randomSolution(): Float32Array {
    const res = new Float32Array(this.inputOpts.dim);

    for (let i = 0; i < this.inputOpts.dim; ++i)
      res[i] = jStat.uniform.sample(this.inputOpts.mins[i], this.inputOpts.maxs[i]);

    return res;
  }

  /** Crossover: blend */
  public crossover(x1: Float32Array, x2: Float32Array): Float32Array {
    const alpha = Math.random();
    const res = new Float32Array(this.inputOpts.dim);

    for (let i = 0; i < this.inputOpts.dim; ++i)
      res[i] = alpha * x1[i] + (1 - alpha) * x2[i];

    return res;
  }

  /** Mutation */
  public mutate(x: Float32Array, sigma: number=0.1): Float32Array {
    const res = new Float32Array(this.inputOpts.dim);

    for (let i = 0; i < this.inputOpts.dim; ++i) {
      res[i] = clip(
        x[i] + jStat.normal.sample(0, sigma),
        this.inputOpts.mins[i],
        this.inputOpts.maxs[i],
      );
    }

    return res;
  }

  /** MOEA/D algorithm */
  public perform(): MoeadOutput[] {
    const nWeights = this.methodOpts.nWeights;
    const weights = this.generateWeights();
    console.log('Weights: ', weights);

    const population = new Array<Float32Array>(nWeights + 1);

    for (let i = 0; i <= nWeights; ++i)
      population[i] = this.randomSolution();

    console.log('Population: ', population);

    const neighborhoods = weights.map((w) =>
      weights
        .map((w2, idx) => ({dist: euclideanDistance(w, w2), idx}))
        .sort((a, b) => a.dist - b.dist)
        .slice(0, this.methodOpts.neighbors)
        .map((e) => e.idx),
    );

    console.log('Neighborhoods: ', neighborhoods);

    // Evolution loop
    for (let gen = 0; gen < this.methodOpts.generations; ++gen) {
      for (let i = 0; i <= nWeights; i++) {
        const nbr = neighborhoods[i];
        const [p1, p2] = pickTwo(nbr);
        let child = this.crossover(population[p1], population[p2]);
        if (Math.random() < this.methodOpts.mutationRate)
          child = this.mutate(child);

        // Update neighborhood
        for (const j of nbr) {
          if (this.weightedSum(child, weights[j]) < this.weightedSum(population[j], weights[j]))
            population[j] = child;
        }
      }
    }

    return population.map((x) => ({
      point: x,
      objective: this.objective(x),
    }));
  }
}; // Moead

const func = (x: Float32Array): Float32Array => {
  return new Float32Array([
    x[0]**2 + x[1]**2 + x[2]**2,
    (x[0] - 2)**2 + (x[1] - 1)**2,
    Math.sqrt((x[0] - 1)**2 + (x[2] - 2)**2),
  ]);
};

const dim = 3;
const inputOpts: InputOptions = {
  dim: dim,
  mins: new Float32Array([0, 0, 0]),
  maxs: new Float32Array([1, 2, 3]),
};
const moead = new Moead(func, inputOpts, dim);
console.log(moead);

// const x = new Float32Array([1, 2]);

// console.log(moead.objective(x));

// const w = moead.generateWeights(2);
// console.log(w);

// console.log(moead.weightedSum(x, w[1]));

// console.log('Random solution:', moead.randomSolution());

// console.log('Crossover:', moead.crossover(
//   new Float32Array([0, 10]),
//   new Float32Array([1, 20]),
// ));

//console.log('Mutation:', moead.mutate(new Float32Array([0.2, 0.3])));

const solution = moead.perform();

for (const sol of solution)
  console.log(sol.point.toString(), ',', sol.objective.toString());

console.log(euclideanDistance(
  new Float32Array([1, 2, 3]),
  new Float32Array([11, 12, 13]),
));
