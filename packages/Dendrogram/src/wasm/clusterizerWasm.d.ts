export declare function getClustersFromDistMatWasm
    (distmat:Float32Array, n: number, method: number): Promise<ClusterMatrix>;
export declare type ClusterMatrix = {
    mergeRow1:Int32Array;
    mergeRow2:Int32Array;
    heightsResult:Float32Array;
}
