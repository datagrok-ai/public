import {tTest, uTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {category, test, expect, expectArray, expectFloat,
  expectExceptionAsync} from '@datagrok-libraries/test/src/test';

// References are computed from the textbook formulae, not read back off this library: sample (n-1)
// variances, Welch-Satterthwaite dof, and t = (m1 - m2) / sqrt(v1/n1 + v2/n2).
const A = [1, 2, 3, 4, 5];
const B = [2, 4, 6, 8, 10];
const EPS = 1e-12;

// R's p.adjust(p, method = 'BH'): the running minimum, taken from the largest p downwards, of
// n/(i+1) * p_i. Independent of alpha - only `reject` depends on it.
const P = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216];
const ADJUSTED = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.10571428, 0.216, 0.216, 0.216];

category('statistics', () => {
  test('tTest: the Welch branch uses the sample variance', async () => {
    const r = tTest(A, B);
    expectFloat(r['p-value'], 0.10753119493044197, EPS);
    expectFloat(r['p-value less'], 0.053765597465220985, EPS);
    expectFloat(r['Mean difference']!, -3, EPS);
    // the population variance, which is what jStat returns without the flag, puts this at 0.0790
    expect(r['p-value'] > 0.1);
  });

  test('tTest: the equal-variance branch divides by the pooled standard deviation', async () => {
    const r = tTest(A, B, false, true);
    expectFloat(r['p-value'], 0.09434977284244564, EPS);
    expectFloat(r['p-value less'], 0.04717488642122282, EPS);
  });

  test('tTest: the known-deviation branch reports a tail probability', async () => {
    const r = tTest(A, B, true);
    expectFloat(r['p-value'], 0.05777957112359733, EPS);
    expectFloat(r['p-value less'], 0.028889785561798664, EPS);
    // a density is bounded by 0.3989 and is blind to the direction; a cdf is neither
    expectFloat(tTest(B, A, true)['p-value less'], 1 - r['p-value less'], EPS);
  });

  test('tTest: the tails complement each other and the two-tailed p is twice the smaller', async () => {
    const results = [tTest(A, B), tTest(A, B, false, true), tTest(A, B, true)];
    for (const r of results) {
      expectFloat(r['p-value less'] + r['p-value more'], 1, EPS);
      expectFloat(r['p-value'], 2 * Math.min(r['p-value less'], r['p-value more']), EPS);
    }
  });

  test('tTest: a sample of fewer than two values is refused', async () => {
    await expectExceptionAsync(async () => void tTest([1], B), (e) => `${e}`.includes('Wrong sample size'));
    await expectExceptionAsync(async () => void tTest(A, []), (e) => `${e}`.includes('Wrong sample size'));
  });

  test('uTest: no ties, hand-computable', async () => {
    // ranks 1..10 with x on the odd ones: R1 = 25, U1 = 10, mu = 12.5, s = sqrt(25/12 * 11)
    const r = uTest([1, 3, 5, 7, 9], [2, 4, 6, 8, 10]);
    expectFloat(r['p-value'], 0.6761033140231467, EPS);
    expectFloat(r['p-value less'], 0.33805165701157336, EPS);
    expectFloat(r['p-value more'], 0.6619483429884266, EPS);
    expectFloat(r['Median difference']!, -1, EPS);
  });

  test('uTest: ties are corrected with the sum of t^3 - t', async () => {
    // xy ranks [1,3,3,5,8,3,6.5,6.5,9,10]: one group of 3 and one of 2, so the term is 24 + 6 = 30.
    // Summing the group sizes instead gives 10 - the sample size - and a p-value of 0.1416.
    const r = uTest([1, 2, 2, 3, 5], [2, 4, 4, 6, 7]);
    expectFloat(r['p-value'], 0.1375638939099033, EPS);
    expectFloat(r['p-value less'], 0.06878194695495166, EPS);
  });

  test('uTest: the one-tailed fields are one-tailed', async () => {
    const r = uTest([1, 3, 5, 7, 9], [2, 4, 6, 8, 10]);
    expect(r['p-value less'] !== r['p-value more']);
    expectFloat(r['p-value less'] + r['p-value more'], 1, EPS);
    expectFloat(r['p-value'], 2 * Math.min(r['p-value less'], r['p-value more']), EPS);

    const swapped = uTest([2, 4, 6, 8, 10], [1, 3, 5, 7, 9]);
    expectFloat(swapped['p-value less'], r['p-value more'], EPS);
    expectFloat(swapped['p-value more'], r['p-value less'], EPS);
    expectFloat(swapped['p-value'], r['p-value'], EPS);
  });

  test('uTest: the two-tailed p never exceeds 1', async () => {
    // the continuity correction used to push U past the mean of its own null distribution
    const r = uTest([1, 2, 3, 4], [1, 2, 3, 4]);
    expectFloat(r['p-value'], 1, EPS);
    expectFloat(r['p-value less'], 0.5, EPS);
  });

  test('uTest: a sample of fewer than two values is refused', async () => {
    await expectExceptionAsync(async () => void uTest([1], [1, 2, 3]), (e) => `${e}`.includes('Wrong sample size'));
    await expectExceptionAsync(async () => void uTest([1, 2, 3], []), (e) => `${e}`.includes('Wrong sample size'));
  });

  test('fdrcorrection: adjusted p-values match Benjamini-Hochberg', async () => {
    const [, adjusted] = fdrcorrection(new Float32Array(P), 0.05, 'i');
    for (let i = 0; i < P.length; i++)
      expectFloat(adjusted[i], ADJUSTED[i], 1e-6);
  });

  test('fdrcorrection: every hypothesis at or below the line is rejected', async () => {
    const alpha = 0.05;
    const [reject, adjusted] = fdrcorrection(new Float32Array(P), alpha, 'i');
    expectArray(reject, [true, true, false, false, false, false, false, false, false, false]);
    // the largest rejected hypothesis used to be dropped, leaving reject out of step with adjusted
    for (let i = 0; i < P.length; i++)
      expect(reject[i], adjusted[i] <= alpha, `hypothesis ${i}`);
  });

  test('fdrcorrection: a lone passing hypothesis is rejected', async () => {
    const [reject] = fdrcorrection(new Float32Array([0.001, 0.9, 0.9, 0.9]), 0.05, 'i');
    expectArray(reject, [true, false, false, false]);
  });

  test('fdrcorrection: rejections come back in the caller order', async () => {
    const shuffled = [0.216, 0.008, 0.212, 0.001, 0.042, 0.041, 0.074, 0.039, 0.205, 0.06];
    const [reject, adjusted] = fdrcorrection(new Float32Array(shuffled), 0.05, 'i');
    expectArray(reject, [false, true, false, true, false, false, false, false, false, false]);
    for (let i = 0; i < shuffled.length; i++)
      expectFloat(adjusted[i], ADJUSTED[P.indexOf(shuffled[i])], 1e-6);
  });

  test('fdrcorrection: nothing is rejected when nothing passes', async () => {
    const [reject] = fdrcorrection(new Float32Array([0.4, 0.5, 0.6]), 0.05, 'i');
    expectArray(reject, [false, false, false]);
  });
});
