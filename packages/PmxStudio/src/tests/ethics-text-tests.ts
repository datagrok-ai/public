import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {EthicsFacts, ethicsText} from '../ethics-text';

const arm = (name: string, n: number | null, power = '80%', effect = 'TGI 60%') =>
  ({name, therapy: 'chemotherapy', effect, n, power});

const base: EthicsFacts = {arms: [arm('Cpd-X 5 mg/kg', 8)], line: 'MC38', papers: 12, target: '80%', alpha: '0.05',
  alphaPer: '0.025', correction: 'holm', nFixed: false, minGroup: 5, preLoss: '10%', postLoss: '5%',
  confidence: '90%', lossRule: 'expected', implant: 13, randomise: 11, nTotal: 8, humane: 2000, diameter: '1.6',
  lastDayChosen: false, end: 21, endpoint: 1500};
const two = [arm('A', 8), arm('B', 8)];

// the approved wording, piece by piece: each test pins the whole text it gives
const BRACKET = '[Explain why an effect of TGI 60% is expected or is the smallest effect of interest.]';
const SIMULATED = 'This is the smallest group size that reaches 80% power, estimated by simulating tumor growth from ' +
  'the controls of 12 published MC38 studies.';
const FEWER = 'Fewer animals per group would fall below 80% power and risk an inconclusive study that may need to be ' +
  'repeated.';
const LOST = 'About 10% of implanted animals are expected to be lost before randomization (no tumor take, tumor ' +
  'outside the volume window) and 5% after it (ulceration, deaths, technical failures).';
const IMPLANT = 'We therefore implant 13, randomize 11 and expect 8 analyzable animals;';
const reduction = (one: string) => `Reduction: ${one}the size is based on the published growth of this tumor model, ` +
  'which avoids using animals in a pilot study for that purpose, and the extra animals only cover the expected ' +
  'losses. An animal reaches its humane endpoint when the tumor exceeds 2000 mm³ (equivalent to a mean diameter of ' +
  'about 1.6 cm); the study ends';
const COMPUTED_END = 'when the median control reaches 1500 mm³, expected around day 21.';
const ONE_ARM = '2 groups: a control and Cpd-X 5 mg/kg (chemotherapy, expected TGI 60%); each treated group is ' +
  'compared with the control.';
const TWO_ARMS = '3 groups: a control and A (chemotherapy, expected TGI 60%), B (chemotherapy, expected TGI 60%); ' +
  'each treated group is compared with the control.';
const request1 = (n: number) => `We request ${n} analyzable animals per treated group and ${n} controls, which ` +
  'gives 80% power to detect TGI 60% at a two-sided significance level of 0.05.';
const REQUEST_2 = 'We request 8 analyzable animals per treated group, with a control group as large as the largest ' +
  'treated group, which gives 80% power to detect TGI 60% at a two-sided significance level of 0.05.';
const paragraph2 = (loss: string, one: string, end: string) =>
  `${LOST} ${IMPLANT} ${loss}. ${reduction(one)} ${end}`;
const AVERAGE = 'the counts cover the average loss';

category('Ethics text', () => {
  test('size computed by simulation, one treated arm', async () => {
    expect(ethicsText(base), `${ONE_ARM} ${request1(8)} ${BRACKET} ${SIMULATED} ${FEWER}\n\n` +
      paragraph2(AVERAGE, '', COMPUTED_END));
  });

  test('size at the minimum of five', async () => {
    expect(ethicsText({...base, arms: [arm('Cpd-X 5 mg/kg', 5)]}), `${ONE_ARM} ${request1(5)} ${BRACKET} ` +
      `${SIMULATED} We use five per group as a minimum, so that losing one animal still leaves a usable ` +
      `comparison.\n\n${paragraph2(AVERAGE, '', COMPUTED_END)}`);
  });

  test('size set by hand', async () => {
    expect(ethicsText({...base, nFixed: true}), `${ONE_ARM} ${request1(8)} ${BRACKET} The group size was chosen by ` +
      `us; the power stated above is what this size provides.\n\n${paragraph2(AVERAGE, '', COMPUTED_END)}`);
  });

  test('size not reachable', async () => {
    expect(ethicsText({...base, arms: [arm('Cpd-X 5 mg/kg', null, '—')], nTotal: null}),
      `${ONE_ARM} No group size up ` +
      `to 40 reaches 80% power for TGI 60%; the study is not adequately powered to detect it. ${BRACKET}\n\n` +
      `${LOST} ${reduction('')} ${COMPUTED_END}`);
  });

  test('several arms, Holm', async () => {
    expect(ethicsText({...base, arms: two}), `${TWO_ARMS} ${REQUEST_2} ${BRACKET} ${SIMULATED} ${FEWER} With 2 ` +
      'comparisons, Holm correction is applied; each group is sized at α = 0.025.\n\n' +
      paragraph2(AVERAGE, 'one control group serves both treated groups; ', COMPUTED_END));
  });

  test('several arms, Bonferroni, loss by confidence, end chosen', async () => {
    expect(ethicsText({...base, arms: two, correction: 'bonferroni', lossRule: 'confidence', lastDayChosen: true}),
      `${TWO_ARMS} ${REQUEST_2} ${BRACKET} ${SIMULATED} ${FEWER} With 2 comparisons, Bonferroni correction is ` +
      'applied; each group is sized at α = 0.025.\n\n' + paragraph2('each count is enough with at least 90% ' +
        'probability', 'one control group serves both treated groups; ', 'on day 21, chosen.'));
  });

  test('several arms, no correction', async () => {
    expect(ethicsText({...base, arms: two, correction: 'none'}), `${TWO_ARMS} ${REQUEST_2} ${BRACKET} ${SIMULATED} ` +
      `${FEWER} Each of the 2 comparisons is tested at α = 0.05, without correction for multiple comparisons.\n\n` +
      paragraph2(AVERAGE, 'one control group serves both treated groups; ', COMPUTED_END));
  });

  test('several arms that differ', async () => {
    expect(ethicsText({...base, arms: [arm('A', 8), arm('B', 12, '82%', 'TGI 45%')]}), `${TWO_ARMS.replace(
      'B (chemotherapy, expected TGI 60%)', 'B (chemotherapy, expected TGI 45%)')} We request 8 analyzable animals ` +
      'for A (80% power to detect TGI 60%), 12 for B (82% power to detect TGI 45%), with a control group as large ' +
      'as the largest treated group, at a two-sided significance level of 0.05. [Explain why effects of TGI 60% and ' +
      `TGI 45% are expected or are the smallest effects of interest.] ${SIMULATED} ${FEWER} With 2 comparisons, ` +
      'Holm correction is applied; each group is sized at α = 0.025.\n\n' +
      paragraph2(AVERAGE, 'one control group serves both treated groups; ', COMPUTED_END));
  });

  test('one arm without a size beside two with one; all three share the control', async () => {
    expect(ethicsText({...base, nTotal: null, arms: [arm('A', 8), arm('B', null, '—', 'TGI 30%'),
      arm('C', 12, '82%', 'TGI 45%')]}), '4 groups: a control and A (chemotherapy, expected TGI 60%), B ' +
      '(chemotherapy, expected TGI 30%), C (chemotherapy, expected TGI 45%); each treated group is compared with the ' +
      'control. We request 8 analyzable animals for A (80% power to detect TGI 60%), 12 for C (82% power to detect ' +
      'TGI 45%), with a control group as large as the largest treated group, at a two-sided significance level of ' +
      '0.05. No group size up to 40 reaches 80% power for TGI 30%; the study is not adequately powered to detect it. ' +
      '[Explain why effects of TGI 60%, TGI 30% and TGI 45% are expected or are the smallest effects of interest.] ' +
      `${SIMULATED} ${FEWER} With 3 comparisons, Holm correction is applied; each group is sized at α = 0.025.\n\n` +
      `${LOST} ${reduction('one control group serves all 3 treated groups; ')} ${COMPUTED_END}`);
  });

  test('two arms without a size; one arm at the minimum beside a larger one', async () => {
    expect(ethicsText({...base, nTotal: null, arms: [arm('A', null, '—'), arm('B', null, '—', 'TGI 45%')]}),
      `${TWO_ARMS.replace('B (chemotherapy, expected TGI 60%)', 'B (chemotherapy, expected TGI 45%)')} No group ` +
      'size up to 40 reaches 80% power for TGI 60% and TGI 45%; the study is not adequately powered to detect them. ' +
      '[Explain why effects of TGI 60% and TGI 45% are expected or are the smallest effects of interest.] With 2 ' +
      'comparisons, Holm correction is applied; each group is sized at α = 0.025.\n\n' +
      `${LOST} ${reduction('one control group serves both treated groups; ')} ${COMPUTED_END}`);
    expect(ethicsText({...base, arms: [arm('A', 5), arm('B', 12, '82%', 'TGI 45%')]}).includes(`${SIMULATED} We use ` +
      'five per group as a minimum, so that losing one animal still leaves a usable comparison. With 2'), true);
  });
});
