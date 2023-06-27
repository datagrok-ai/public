import {mmDistanceFunctionType} from './types';

export function hamming(): mmDistanceFunctionType {
  return (seq1: string, seq2: string) => {
    // hamming distance should only be used with same size strings,
    // but still, lets add a check and if they are not same length add the difference to the result
    let diff = 0;
    if (seq1.length !== seq2.length)
      diff = Math.abs(seq1.length - seq2.length);

    let result = 0;
    for (let i = 0; i < Math.min(seq1.length, seq2.length); i++) {
      if (seq1[i] !== seq2[i])
        result++;
    }
    result += diff;
    result /= Math.max(seq1.length, seq2.length);
    return result;
  };
}
