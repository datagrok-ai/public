/* eslint-disable valid-jsdoc */
// The MOEA/D method for multi-objective optimization: https://ieeexplore.ieee.org/document/4358754

//@ts-ignore: no types
import * as jStat from 'jstat';
import {DEFAULT_SETTINGS, Func, InputOptions, MoeadOptions, SLEEP_TIME} from './defs';
import {clip, euclideanDistance, pickTwo} from './utils';

function sleep(ms: number) {
  return new Promise((resolve, reject) => setTimeout(resolve, ms));
};

/** The MOEA/D multi-objective optimizer */
export class Moead {
  private objective: Func;
  private inputOpts: InputOptions;
  private outputDim: number;
  private methodOpts: MoeadOptions;
  private progrIndicator: any = undefined;

  constructor(
    objective: Func,
    inputOpts: InputOptions,
    outputDim: number,
    methodOpts: MoeadOptions = DEFAULT_SETTINGS,
    progrIndicator?: any,
  ) {
    this.objective = objective;
    this.inputOpts = inputOpts;
    this.outputDim = outputDim;
    this.methodOpts = methodOpts;
    this.progrIndicator = progrIndicator;
  }

  /** Scalarization: weighted sum */
  private async weightedSum(x: Float32Array, w: Float32Array): Promise<number> {
    let sum = 0;
    const out = await this.objective(x);

    for (let i = 0; i < this.outputDim; ++i)
      sum += w[i] * out[i];

    return sum;
  }

  /** Generate uniform weight vectors */
  private generateWeights(): Float32Array[] {
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
  private randomSolution(): Float32Array {
    const res = new Float32Array(this.inputOpts.dim);

    for (let i = 0; i < this.inputOpts.dim; ++i)
      res[i] = jStat.uniform.sample(this.inputOpts.mins[i], this.inputOpts.maxs[i]);

    return res;
  }

  /** Crossover: blend */
  private crossover(x1: Float32Array, x2: Float32Array): Float32Array {
    const alpha = Math.random();
    const res = new Float32Array(this.inputOpts.dim);

    for (let i = 0; i < this.inputOpts.dim; ++i)
      res[i] = alpha * x1[i] + (1 - alpha) * x2[i];

    return res;
  }

  /** Mutation */
  private mutate(x: Float32Array, sigma: number=0.1): Float32Array {
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
  public async perform(): Promise<Float32Array[]> {
    const nWeights = this.methodOpts.nWeights;
    const weights = this.generateWeights();
    const population = new Array<Float32Array>(nWeights + 1);

    for (let i = 0; i <= nWeights; ++i)
      population[i] = this.randomSolution();

    const neighborhoods = weights.map((w) =>
      weights
        .map((w2, idx) => ({dist: euclideanDistance(w, w2), idx}))
        .sort((a, b) => a.dist - b.dist)
        .slice(0, this.methodOpts.neighbors)
        .map((e) => e.idx),
    );

    const iter = this.methodOpts.generations;
    let percentage = 0;

    // Evolution loop
    for (let gen = 0; gen < iter; ++gen) {
      for (let i = 0; i <= nWeights; i++) {
        const nbr = neighborhoods[i];
        const [p1, p2] = pickTwo(nbr);
        let child = this.crossover(population[p1], population[p2]);
        if (Math.random() < this.methodOpts.mutationRate)
          child = this.mutate(child);

        // Update neighborhood
        for (const j of nbr) {
          if (await this.weightedSum(child, weights[j]) < await this.weightedSum(population[j], weights[j]))
            population[j] = child;
        }
      }

      if (this.progrIndicator !== undefined) {
        if (this.progrIndicator.canceled)
          throw new Error('Optimization interrupted by the user');

        percentage = Math.floor(100 * gen / iter);
        this.progrIndicator.update(percentage, `Optimization... (${percentage}%)`);
        await sleep(SLEEP_TIME);
      }
    }

    return population;
  }
}; // Moead
