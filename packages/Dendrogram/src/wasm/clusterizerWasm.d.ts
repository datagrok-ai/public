import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';

export declare function getClustersFromDistMatWasm
    (distmat:Float32Array, n: number, method: number): Promise<ClusterMatrix>;
