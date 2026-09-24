import {sampleReactionsWasm} from './sampler-wrapper';
import type {PrecomputedWarmup} from './sampler-wrapper';

onmessage = async (message) => {
  const {model, samples, thinning, precomputedWarmup} = message.data as
    {model: any; samples: number; thinning: number; precomputedWarmup?: PrecomputedWarmup};
  try {
    const modelResult = await sampleReactionsWasm(model, samples, thinning, precomputedWarmup);
    postMessage(modelResult);
  } catch (e: any) {
    postMessage({error: e?.message ?? String(e)});
  }
};
