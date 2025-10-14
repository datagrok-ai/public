import { WorkerCobraSolver } from '.';
import type { CobraModelData } from '../../escher_src/src/ts/types';
import {exportCppSampler} from './sampler.js';

class SamplerWasm {
    private static _wasmInstance: any | null = null;
    static async getInstance() {
        if (this._wasmInstance == null) {
            const wasmUrl = new URL('./sampler.wasm', import.meta.url).href;
            const wasmPath = wasmUrl.substring(0, wasmUrl.lastIndexOf('/') + 1) + 'sampler.wasm';
            try {
                this._wasmInstance = await exportCppSampler();
            } catch (e) {
                try {
                    this._wasmInstance = await exportCppSampler({locateFile: () => 'sampler.wasm'});
                }
                catch (e) {
                    console.error(e);
                    throw new Error('Unable to load wasm file for dbscan');
                }
            }
        }
        return this._wasmInstance;
    }
}

export async function sampleReactionsWasm(mp: CobraModelData, samplesCount: number = 1000, thinning: number = 20): Promise<Float32Array> {
    
    const wasmInstance = await SamplerWasm.getInstance();
    const getSamples = wasmInstance.cwrap('sample', 'number', ['number', 'number','number','number','number','number','number','number','number','number']);
    // arguments:
    // int sampleCount
    // int thinning
    // int reactionCount
    // int metaboliteCount
    // Float32Array lbs
    // Float32Array ubs
    // Float32Array sData (S matrix or the stoichiometric matrix. each row represents a metabolite, each column a reaction)
    // int initOptRows - number of rows in the initial optimization problem
    // Float32Array initOptData - each row is one sample, each column is a reaction, stacked row-wise
    // Float32Array - pointer to the output array, should be of size sampleCount * reactionCount
    const reactionCount = mp.reactions.length;
    const metaboliteCount = mp.metabolites.length;
    const lbs = new Float32Array(reactionCount).fill(0).map((_, i) => mp.reactions[i].lower_bound ?? 0);
    const ubs = new Float32Array(reactionCount).fill(0).map((_, i) => mp.reactions[i].upper_bound ?? 1000);
    const sData = new Float32Array(metaboliteCount * reactionCount);
    const metaboliteLookup: Map<string, number> = new Map();
    mp.reactions.forEach((r, j) => {
        r.metabolites && Object.entries(r.metabolites).forEach(([metName, coeff]) => {
            let metIndex = metaboliteLookup.get(metName);
            if (metIndex == undefined) {
                metIndex = metaboliteLookup.size;
                metaboliteLookup.set(metName, metIndex);
            }
            sData[metIndex * reactionCount + j] = coeff;
        });
    });

    // run initial optimization to get a starting point
    const optiExtremes = await WorkerCobraSolver.get_extreme_points(mp);
    console.log('got extremes');
    const initOpt = new Float32Array(reactionCount * optiExtremes.solutions.length);
    for (let i = 0; i < optiExtremes.solutions.length; i++)
        initOpt.set(optiExtremes.solutions[i], i * reactionCount);

    const lbsPointer = wasmInstance._malloc(lbs.length * lbs.BYTES_PER_ELEMENT);
    const ubsPointer = wasmInstance._malloc(ubs.length * ubs.BYTES_PER_ELEMENT);
    const sDataPointer = wasmInstance._malloc(sData.length * sData.BYTES_PER_ELEMENT);
    const initOptPointer = wasmInstance._malloc(initOpt.length * initOpt.BYTES_PER_ELEMENT);
    const outputPointer = wasmInstance._malloc(samplesCount * reactionCount * Float32Array.BYTES_PER_ELEMENT * 2);
    
    const lbsHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, lbsPointer, lbs.length * lbs.BYTES_PER_ELEMENT);
    const ubsHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, ubsPointer, ubs.length * ubs.BYTES_PER_ELEMENT);
    const sDataHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, sDataPointer, sData.length * sData.BYTES_PER_ELEMENT);
    const initOptHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, initOptPointer, initOpt.length * initOpt.BYTES_PER_ELEMENT);
    const outputHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, outputPointer, samplesCount * reactionCount * Float32Array.BYTES_PER_ELEMENT);
    lbsHeap.set(new Uint8Array(lbs.buffer));
    ubsHeap.set(new Uint8Array(ubs.buffer));
    sDataHeap.set(new Uint8Array(sData.buffer));
    initOptHeap.set(new Uint8Array(initOpt.buffer));
    // outputHeap.fill(0);

    getSamples(samplesCount, thinning, reactionCount, metaboliteCount,
        lbsPointer, ubsPointer, sDataPointer, optiExtremes.solutions.length, initOptPointer, outputPointer
    )

    const getSamplesResult = new Float32Array(outputHeap.buffer, outputHeap.byteOffset, samplesCount * reactionCount);

    wasmInstance._free(lbsPointer);
    wasmInstance._free(ubsPointer);
    wasmInstance._free(sDataPointer);
    wasmInstance._free(initOptPointer);
    wasmInstance._free(outputPointer);
    
    return getSamplesResult;
}