/* eslint-disable valid-jsdoc */
// The MOEA/D method for multi-objective optimization: https://ieeexplore.ieee.org/document/4358754

/**
 import numpy as np
import matplotlib.pyplot as plt

# Objective functions (now depend on x1 and x2)
def f1(x):
    return x[0]**2 + x[1]**2

def f2(x):
    return (x[0] - 2)**2 + (x[1] - 1)**2

# Scalarization: weighted sum
def weighted_sum(x, w):
    return w[0] * f1(x) + w[1] * f2(x)

# Generate uniform weight vectors
def generate_weights(n_weights):
    weights = []
    for i in range(n_weights + 1):
        w1 = i / n_weights
        w2 = 1 - w1
        weights.append((w1, w2))
    return np.array(weights)

# Initialize a random 2D solution in [0, 2] × [0, 2]
def random_solution():
    return np.random.uniform(0, 2, size=2)

# Crossover: blend
def crossover(x1, x2):
    alpha = np.random.rand()
    return alpha * x1 + (1 - alpha) * x2

# Mutation
def mutate(x, sigma=0.1):
    return np.clip(x + np.random.normal(0, sigma, size=x.shape), 0, 2)

# MOEA/D algorithm
def moead(n_weights=100, generations=100, neighbors=10, mutation_rate=0.2):
    weights = generate_weights(n_weights)
    population = np.array([random_solution() for _ in range(n_weights + 1)])

    # Precompute neighborhoods
    distances = np.linalg.norm(weights[:, None, :] - weights[None, :, :], axis=2)
    neighborhoods = np.argsort(distances, axis=1)[:, :neighbors]

    for gen in range(generations):
        for i in range(n_weights + 1):
            nbr = neighborhoods[i]
            p1, p2 = np.random.choice(nbr, 2, replace=False)
            x1, x2 = population[p1], population[p2]

            # Create offspring
            offspring = crossover(x1, x2)
            if np.random.rand() < mutation_rate:
                offspring = mutate(offspring)

            # Update neighborhood
            for j in nbr:
                if weighted_sum(offspring, weights[j]) < weighted_sum(population[j], weights[j]):
                    population[j] = offspring

    return population, weights

# Run MOEA/D
solutions, weights = moead()

# Evaluate objectives
f1_vals = np.array([f1(x) for x in solutions])
f2_vals = np.array([f2(x) for x in solutions])

# Plot Pareto front
plt.figure(figsize=(7, 5))
plt.scatter(f1_vals, f2_vals, c='blue', label='MOEA/D Approximation')
plt.xlabel('f1(x)')
plt.ylabel('f2(x)')
plt.title('Pareto Front for 2D Variables')
plt.grid(True)
plt.legend()
plt.show()

 */

//@ts-ignore: no types
import * as jStat from 'jstat';
import {DEFAULT_SETTINGS, Func, InputOptions, MoeadOptions} from './defs';
import {clip, euclideanDistance} from './utils';

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
  public perform() {
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

    //np.array([random_solution() for _ in range(n_weights + 1)])
  }
}; // Moead

const func = (x: Float32Array): Float32Array => {
  return new Float32Array([
    x[0]**2 + x[1]**2,
    (x[0] - 2)**2 + (x[1] - 1)**2,
  ]);
};

const dim = 2;
const inputOpts: InputOptions = {
  dim: dim,
  mins: new Float32Array([0, 0]),
  maxs: new Float32Array([2, 2]),
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

console.log(moead.perform());

console.log(euclideanDistance(
  new Float32Array([1, 2, 3]),
  new Float32Array([11, 12, 13]),
));
