import {mmDistanceFunctionArgs, mmDistanceFunctionType} from './types';

export function hamming(args: Partial<mmDistanceFunctionArgs> = {}): mmDistanceFunctionType {
  function getDistanceF(): (a: string, b: string) => number {
    if (!args || !args.scoringMatrix || !args.alphabetIndexes)
      return (a: string, b: string) => a === b ? 0 : 1;
    if (args.scoringMatrix.length !== Object.keys(args.alphabetIndexes).length)
      throw new Error('Scoring matrix and alphabet indexes should have the same length');
    const indexes = args.alphabetIndexes;
    const matrix = args.scoringMatrix;
    //const matrixMap = new Map<string, Map<string, number>>();
    //const map2: any = {};
    const minCharCode = Math.min(...Object.keys(indexes).map((k) => k.charCodeAt(0))) + 1;

    const scorringArray = new Float32Array((matrix.length + minCharCode) * (matrix.length + minCharCode));
    Object.entries(indexes).forEach(([key, index]) => {
      //matrixMap.set(key, new Map<string, number>());
      //map2[key] = {};
      const matrixRow = matrix[index];
      Object.entries(indexes).forEach(([key2, index2]) => {
        //matrixMap.get(key)!.set(key2, matrixRow[index2]);
        scorringArray[key.charCodeAt(0) * matrix.length + key2.charCodeAt(0)] = matrixRow[index2];
        //map2[key][key2] = matrixRow[index2];
      });
    });
    return (a: string, b: string) => {
      return scorringArray[a.charCodeAt(0) * matrix.length + b.charCodeAt(0)];
    };
  }
  const distanceF = getDistanceF();

  const threshold = args?.threshold ?? 0;

  return (seq1: string, seq2: string) => {
    // hamming distance should only be used with same size strings,
    // but still, lets add a check and if they are not same length add the difference to the result
    let diff = 0;
    const s1l = seq1.length;
    const s2l = seq2.length;
    const thresholdLimit = Math.max(s1l, s2l) * (1 - threshold);
    if (s1l !== s2l)
      diff = Math.abs(s1l - s2l);

    let result = 0;
    for (let i = 0; i < Math.min(s1l, s2l); i++) {
      if (seq1[i] !== seq2[i]) {
        result += distanceF(seq1[i], seq2[i]);
        if (result > thresholdLimit)
          return 1;
      }
    }
    result += diff;
    result /= Math.max(s1l, s2l);
    return result;
  };
}
