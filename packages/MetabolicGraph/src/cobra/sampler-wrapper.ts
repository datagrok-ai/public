import {exportCppSampler} from './sampler.js';
import type {CobraModelData} from '../../escher_src/src/ts/types';

/**
 * cobra OptGP warmup points computed elsewhere (the ComputeExtremePoints Python script): each point
 * lists the forward variable of every reaction, then every reverse variable, in `reactionNames` order.
 */
export type PrecomputedWarmup = {
  reactionNames: string[];
  points: number[][];
};

/** The seed of the Python route (OptGpSampling.py), so both routes return the same samples. */
export const SAMPLER_SEED = 42;

// negative return codes of sample() / sample_with_warmup() in sampler.cpp
const SAMPLER_ERRORS: {[code: number]: string} = {
  [-1]: 'the feasible flux space is a single point, there is nothing to sample',
  [-2]: 'cannot sample an inhomogeneous problem with only 2 search directions',
  [-3]: 'the model seems numerically unstable, sampling could not proceed',
  [-4]: 'infeasible reaction bounds: no flux distribution satisfies them',
  [-5]: 'could not set up the warmup problem',
};

// Referencing the binary makes webpack emit sampler.wasm next to the bundles, where the module finds it.
new URL('./sampler.wasm', import.meta.url);

class SamplerWasm {
  private static _instance: Promise<any> | null = null;
  static getInstance(): Promise<any> {
    this._instance ??= exportCppSampler().catch(() => exportCppSampler({locateFile: () => 'sampler.wasm'}))
      .catch((e: any) => {
        this._instance = null;
        console.error(e);
        throw new Error('Unable to load the sampler WebAssembly module');
      });
    return this._instance!;
  }
}

/** Warmup points reordered to the model's reactions, flattened for the sampler. */
function warmupInModelOrder(mp: CobraModelData, warmup: PrecomputedWarmup): Float64Array {
  const n = mp.reactions.length;
  const col = new Map(warmup.reactionNames.map((id, i) => [id, i]));
  const order = mp.reactions.map((r) => {
    const i = col.get(r.id);
    if (i === undefined)
      throw new Error(`Warmup points do not cover reaction ${r.id}`);
    return i;
  });
  const m = warmup.reactionNames.length;
  const out = new Float64Array(warmup.points.length * 2 * n);
  warmup.points.forEach((p, k) => {
    for (let j = 0; j < n; j++) {
      out[k * 2 * n + j] = p[order[j]];
      out[k * 2 * n + n + j] = p[m + order[j]];
    }
  });
  return out;
}

/**
 * Samples flux vectors (stacked row-wise, one row per sample) exactly like cobra's
 * OptGPSampler(model, thinning, processes=1, seed).sample(samplesCount): the WebAssembly module
 * computes cobra's warmup points itself (or takes `precomputedWarmup`) and runs the same chain.
 */
export async function sampleReactionsWasm(mp: CobraModelData, samplesCount: number = 1000, thinning: number = 20,
  precomputedWarmup?: PrecomputedWarmup, seed: number = SAMPLER_SEED): Promise<Float32Array> {
  const wasm = await SamplerWasm.getInstance();
  const n = mp.reactions.length;
  const lbs = Float64Array.from(mp.reactions, (r) => r.lower_bound ?? 0);
  const ubs = Float64Array.from(mp.reactions, (r) => r.upper_bound ?? 1000);
  const inverted = mp.reactions.find((_, i) => !(lbs[i] <= ubs[i]));
  if (inverted)
    throw new Error(`Reaction ${inverted.id} has its lower bound above its upper bound`);

  // one stoichiometric row per metabolite, in model order like cobra's solver problem: the warmup
  // LPs pivot differently otherwise. Metabolites missing from the list are appended.
  const rowOf = new Map<string, number>(mp.metabolites.map((m, i) => [m.id, i]));
  for (const r of mp.reactions) {
    for (const id of Object.keys(r.metabolites ?? {})) {
      if (!rowOf.has(id))
        rowOf.set(id, rowOf.size);
    }
  }
  const m = rowOf.size;
  const S = new Float64Array(m * n);
  mp.reactions.forEach((r, j) => {
    for (const [id, coef] of Object.entries(r.metabolites ?? {}))
      S[rowOf.get(id)! * n + j] = coef;
  });
  const warmup = precomputedWarmup ? warmupInModelOrder(mp, precomputedWarmup) : null;

  // allocate everything first: an allocation can grow the heap and detach earlier views
  const ptrs = {
    lbs: wasm._malloc(lbs.byteLength), ubs: wasm._malloc(ubs.byteLength), S: wasm._malloc(S.byteLength),
    warmup: warmup ? wasm._malloc(warmup.byteLength) : 0,
    result: wasm._malloc(samplesCount * n * Float32Array.BYTES_PER_ELEMENT),
    stats: wasm._malloc(3 * Int32Array.BYTES_PER_ELEMENT),
  };
  try {
    wasm.HEAPF64.set(lbs, ptrs.lbs / 8);
    wasm.HEAPF64.set(ubs, ptrs.ubs / 8);
    wasm.HEAPF64.set(S, ptrs.S / 8);
    if (warmup)
      wasm.HEAPF64.set(warmup, ptrs.warmup / 8);
    const code: number = warmup ?
      wasm._sample_with_warmup(samplesCount, thinning, seed, n, m, ptrs.lbs, ptrs.ubs, ptrs.S,
        warmup.length / (2 * n), ptrs.warmup, ptrs.result, ptrs.stats) :
      wasm._sample(samplesCount, thinning, seed, n, m, ptrs.lbs, ptrs.ubs, ptrs.S, ptrs.result, ptrs.stats);
    if (code < 0)
      throw new Error(`Sampling failed: ${SAMPLER_ERRORS[code] ?? `error ${code}`}`);
    return new Float32Array(wasm.HEAPF32.buffer, ptrs.result, samplesCount * n).slice();
  } finally {
    for (const p of Object.values(ptrs)) {
      if (p)
        wasm._free(p);
    }
  }
}
