import {pass1, pass2, pass3, computeMedian} from '../calculators';
import {
  constantData, rampData, alternatingData, spikeData,
} from './helpers';

describe('pass1', () => {
  it('constant series', () => {
    const p = pass1(constantData, 10);
    expect(p.mean).toBe(5);
    expect(p.min).toBe(5);
    expect(p.max).toBe(5);
    expect(p.std).toBe(0);
    expect(p.variance).toBe(0);
    expect(p.range).toBe(0);
    expect(p.minCount).toBe(10);
    expect(p.maxCount).toBe(10);
  });

  it('ramp', () => {
    const p = pass1(rampData, 10);
    expect(p.min).toBe(0);
    expect(p.max).toBe(9);
    expect(p.firstMinIdx).toBe(0);
    expect(p.firstMaxIdx).toBe(9);
    expect(p.lastMinIdx).toBe(0);
    expect(p.lastMaxIdx).toBe(9);
    expect(p.mean).toBeCloseTo(4.5, 10);
    expect(p.sum).toBe(45);
    expect(p.sumSq).toBe(285);
  });

  it('spike — min/max tracking', () => {
    const p = pass1(spikeData, 10);
    expect(p.min).toBe(0);
    expect(p.max).toBe(10);
    expect(p.firstMinIdx).toBe(0);
    expect(p.firstMaxIdx).toBe(4);
    expect(p.lastMinIdx).toBe(9);
    expect(p.lastMaxIdx).toBe(4);
    expect(p.minCount).toBe(9);
    expect(p.maxCount).toBe(1);
  });
});

describe('pass2', () => {
  it('constant series — all diffs zero, all duplicates', () => {
    const p = pass2(constantData, 10, 5, 0);
    expect(p.sumAbsDiff).toBe(0);
    expect(p.sumSqDiff).toBe(0);
    expect(p.hasDuplicate).toBe(true);
    expect(p.uniqueCount).toBe(1);
    expect(p.reoccurringDatapointCount).toBe(10);
  });

  it('alternating — high diffs, crossings', () => {
    const p = pass2(alternatingData, 10, 0, 0);
    expect(p.crossingCount).toBe(9);
    expect(p.sumAbsDiff).toBe(18);
  });

  it('ramp — unique values, one crossing of 0', () => {
    const p = pass2(rampData, 10, 4.5, 0);
    expect(p.uniqueCount).toBe(10);
    expect(p.hasDuplicate).toBe(false);
    expect(p.crossingCount).toBe(1);
    expect(p.longestStrikeAbove).toBe(5);
    expect(p.longestStrikeBelow).toBe(5);
  });
});

describe('pass3', () => {
  it('ramp — count above/below mean', () => {
    const p = pass3(rampData, 10, 4.5, 2.8722813232690143, [1, 1.5, 2]);
    expect(p.countAboveMean).toBe(5);
    expect(p.countBelowMean).toBe(5);
  });

  it('constant series — nothing above or below', () => {
    const p = pass3(constantData, 10, 5, 0, [1, 1.5, 2]);
    expect(p.countAboveMean).toBe(0);
    expect(p.countBelowMean).toBe(0);
    expect(p.beyondRSigmaCounts).toEqual([0, 0, 0]);
  });

  it('spike — ratio_beyond_r_sigma', () => {
    const p = pass3(spikeData, 10, 1, 3, [1, 1.5, 2]);
    // spike at 10 is 9 away from mean=1, that's 3*std. So > 1*std, > 1.5*std, > 2*std
    expect(p.beyondRSigmaCounts[0]).toBe(1); // > 1*3 = 3
    expect(p.beyondRSigmaCounts[1]).toBe(1); // > 1.5*3 = 4.5
    expect(p.beyondRSigmaCounts[2]).toBe(1); // > 2*3 = 6
  });
});

describe('computeMedian', () => {
  it('even n', () => {
    const buf = new Float64Array(4);
    expect(computeMedian(new Float64Array([1, 2, 3, 4]), 4, buf)).toBe(2.5);
  });

  it('odd n', () => {
    const buf = new Float64Array(5);
    expect(computeMedian(new Float64Array([1, 2, 3, 4, 5]), 5, buf)).toBe(3);
  });

  it('does not mutate input', () => {
    const input = new Int32Array([3, 1, 2]);
    const buf = new Float64Array(3);
    const result = computeMedian(input, 3, buf);
    expect(result).toBe(2);
    expect(input[0]).toBe(3);
    expect(input[1]).toBe(1);
    expect(input[2]).toBe(2);
  });

  it('unsorted data', () => {
    const buf = new Float64Array(5);
    expect(computeMedian(new Float64Array([5, 1, 4, 2, 3]), 5, buf)).toBe(3);
  });
});
