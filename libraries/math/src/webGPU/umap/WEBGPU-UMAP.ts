/* eslint-disable max-len */

import {performMatrixOps} from './fuzzy-simplical-set';
import {computeMembershipStrengths} from './membership-strengths';
import {optimizationLoop} from './optimization-loop';
import {initializeSimplicialSetEmbedding} from './simplical-set-embedding';
import {smoothKNNDistance} from './smooth-knn-distance';

// TODO: FIT ASYNC METHOD WITH CALLBACKS on each step
// TODO: random function passing in params


export class UmapParams {
  nComponents: number = 2;
  gamma: number = 1;
  alpha: number = 1;
  entryLen: number = 0;
  nNeighbours: number = 15;
  spread: number = 1;
  minDist: number = 0.1;
  negativeSampleRate: number = 5;
  localConnectivity: number = 1.0;
  setOpMixRatio: number = 1.0;
}

export class WebGPUUMAP {
  params: UmapParams = new UmapParams();
  knnIndexes: Int32Array[] | Uint32Array[] | number[][] | null = null;
  knnDistances: Float32Array[] | number[][] | null = null;
  rhos: Float32Array | null = null;
  sigmas: Float32Array | null = null;
  membershipStrengths: Float32Array[] | null = null;


  constructor(entryLen: number, params: Partial<UmapParams>) {
    this.params.entryLen = entryLen;
    Object.assign(this.params, params);
  }

  setPrecomputedKNN(knnIndexes: Int32Array[] | Uint32Array[] | number[][], knnDistances: Float32Array[] | number[][]) {
    this.knnIndexes = knnIndexes;
    this.knnDistances = knnDistances;
  }

  // for now, not sure if we should include this, maybe webGPU UMAP should be stand alone.
  // eslint-disable-next-line max-len
  //   async calcKNN(entryList: SupportedEntryTypes[][], distanceMetrics: WEBGPUDISTANCE[], aggregationFunction: WEBGSLAGGREGATION, weights: number[], options: {
  //         [key: string]: any;
  //     }[]) {
  //     try {
  //       const knnRes = await multiColWebGPUKNN(entryList, this.params.nNeighbours, distanceMetrics, aggregationFunction, weights, options);
  //       if (knnRes) {
  //         this.knnIndexes = knnRes.knnIndexes;
  //         this.knnDistances = knnRes.knnDistances;
  //       }
  //     } catch (e) {
  //       console.error(e);
  //     }
  //     if (!this.knnIndexes || !this.knnDistances)
  //       throw new Error('failed to compute knn indexes and distances');
  //   }

  //   async calcKNNSingle(entryList: SupportedEntryTypes[], distanceMetrics: WEBGPUDISTANCE, options?: {
  //         [key: string]: any;
  //     }) {
  //     if (!options)
  //       options = {};

  //     await this.calcKNN([entryList], [distanceMetrics], WEBGSLAGGREGATION.MANHATTAN, [1], [options]);
  //   }


  private async initializeFit() {
    await this.calcSmoothKNNDistances();
    await this.calcMembershipStrengths();
  }

  public async fit() {
    await this.initializeFit();
    const umapParams = await this.initializeSimplicialSetEmbedding();
    if (!umapParams)
      throw new Error('failed to compute umap simplistical set embeddings');

    const {
      head,
      tail,
      epochsPerSample,
      epochsPerNegativeSample,
      a,
      b,
      gamma,
      initialAlpha,
      nComponents,
      nEpochs,
      entryLen
    } = umapParams;

    const resEmbeddings = await optimizationLoop(
      head, tail, entryLen, epochsPerSample, epochsPerNegativeSample, initialAlpha, gamma, a, b, nComponents, nEpochs, entryLen
    );
    return resEmbeddings;
  }

  private async calcSmoothKNNDistances() {
    if (!this.knnIndexes || !this.knnDistances)
      throw new Error('knn indexes and distances must be set before calling fit');
    const res = await smoothKNNDistance(this.knnDistances, this.params.nNeighbours, this.params.localConnectivity, 64);

    if (!res)
      throw new Error('failed to compute smooth knn distances');
    this.rhos = res.resultRhos;
    this.sigmas = res.resultSigmas;
  }

  private async calcMembershipStrengths() {
    if (!this.knnIndexes || !this.knnDistances || !this.rhos || !this.sigmas)
      throw new Error('knn indexes, distances, rhos, and sigmas must be set before calling fit');

    const res = await computeMembershipStrengths(this.knnDistances, this.sigmas, this.rhos);
    if (!res)
      throw new Error('failed to compute membership strengths');
    this.membershipStrengths = res;
  }

  private async fuzzySimplicialSet() {
    if (!this.knnIndexes || !this.knnDistances || !this.rhos || !this.sigmas || !this.membershipStrengths)
      throw new Error('knn indexes, distances, rhos, sigmas, and membership strengths must be set before calling fit');

    const res = await performMatrixOps(this.knnIndexes, this.membershipStrengths, this.params.setOpMixRatio);
    if (!res)
      throw new Error('failed to compute pairwise multiply with transpose');
    return res;
  }

  private async initializeSimplicialSetEmbedding() {
    const fuzzySimplicialSetRes = await this.fuzzySimplicialSet();
    if (!fuzzySimplicialSetRes || !fuzzySimplicialSetRes.res || !fuzzySimplicialSetRes.unionMatrixOffsets || !fuzzySimplicialSetRes.unionSizes)
      throw new Error('failed to compute fuzzy simplicial set');

    const eRes = await initializeSimplicialSetEmbedding(fuzzySimplicialSetRes, this.params.entryLen, this.params.spread, this.params.minDist, this.params.negativeSampleRate);
    if (!eRes)
      throw new Error('failed to compute umap simplistical set embeddings');
    return eRes;
  }
}
