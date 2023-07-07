import {ClusterMatrix} from '../types';

export declare function getClustersFromDistMatWasm
    (distmat:Float32Array, n: number, method: number): Promise<ClusterMatrix>;
