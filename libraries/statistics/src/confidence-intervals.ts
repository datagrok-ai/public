//@ts-ignore: no types
import * as jStat from 'jstat';

export type Confidence = {
  central: number,
  top: number,
  bottom: number
}

/**
 * @param {number[]} x observations vector.
 * @param {number} confidenceLevel 0.05 means that 95% of observations enter the interval.
 * @param {boolean} parametric symmetric parametric interval if true, sample quintiles if false.
 * @return {Confidence} Object containing central tendency estimate and top and bottom of ht interval.
 */
export function getConfidence(x: number[], confidenceLevel: number = 0.05, parametric: boolean = true):
  Confidence {
  if (x.length <= 1) 
    throw "vector is too short";

  if (confidenceLevel >=1 || confidenceLevel <=0)
    throw "incorrect confidence level";

  if (parametric) {
    const average = jStat.mean(x);
    //true flag is obligatory for sample variance
    const sigma = jStat.stdev(x, true);
    //degrees of freedom as n-1 for student distribution
    const interval = jStat.studentt.inv(1 - confidenceLevel/2, x.length - 1)*sigma/Math.sqrt(x.length);

    const res:Confidence = {
      central: average,
      top: average + interval,
      bottom: average - interval
    };  

    return res;
  } else {
    const res:Confidence = {
      central: jStat.median(x),
      top: jStat.percentile(x, 1 - confidenceLevel/2),
      bottom: jStat.percentile(x, confidenceLevel/2)
    };  

    return res;
  }
}