/* eslint-disable camelcase */
/* eslint-disable max-len */
export const DEFAULT_EQUATIONS = `#name: UDel
#tags: model
#description: Bioreactor simulation, ver. Nov 19, 2024

#equations:
  d(X_T)/dt  =   mu * X_V - K_lysis * (X_T - X_V)
  d(X_V)/dt  =  (mu - k_d) * X_V
  d(Glc)/dt  =   dU3
  d(Gln)/dt  =   dU4
  d(Lac)/dt  =   dU5
  d(Amm)/dt  =  - Y_Amm_Gln * dU4
  d(MAb)/dt  =   dMab
  d(ATP)/dt  =   dU7  
  d(Vol)/dt  =   0

#expressions:
  mu_Glc        = Glc / (K_Glc + Glc)
  mu_Gln        = Gln / (K_Gln + Gln)
  Glc_per_Cell  = Glc / X_V
  sigma_Glc     = 1 / (1 + exp(sigm_factor * (Glc_per_Cell - Glc_Cell_Th)))
  mu_Lac        = mu_1_max * sigma_Glc * K_Lac_i / (K_Lac_i + Lac) + mu_2_max * (1 - sigma_Glc) * Lac / (K_Lac + Lac)
  mu_Amm        = K_Amm / (K_Amm + Amm)
  mu            = mu_Gln * mu_Amm * mu_Lac * mu_Glc
  k_d           = k_d_max * k_mu / (mu + k_mu)
  dU3           = - mu / Y_X_Glc * X_V - m_Glc * X_V
  dU4           = - mu / Y_X_Gln * X_V - m_Gln * X_V
  dU5           = - Y_Lac_Glc * dU3 * sigma_Glc - mu / Y_X_Lac * X_V * (1 - sigma_Glc)
  dU7           = - sigma_Glc * 30 * dU3 - (1 - sigma_Glc) * 12 * dU5 - 15 * dU4
  dMab          = (1/1342 * dU7 - 1/66 * dU3 - 1/62 * dU4)
  Prod          = dMab
  denom         = X_V * W_dry_scale * W_dry
  SPr           = dMab / denom
  SSUR          = -dU3 / denom

#argument: t
  initial = 0     {units: h; category: Time; min: 0; max: 10}   [Begin of simulation interval]
  final = 250   {units: h; category: Time; min: 20; max: 750} [Feeding interval interval]
  step = 3      {units: h; category: Time; min: 0.5; max: 24} [Time step of simulation]

#output:
  t
  X_T
  X_V
  Glc
  Gln
  Lac
  Amm
  MAb
  ATP
  Vol
  Prod
  SPr
  SSUR

#inits:
  X_T = 1      {category: Initials; units: mln. cells / L}   [Initial total cell density]
  X_V = 1      {category: Initials; units: mln. cells / L}   [Initial viable cell density]
  Glc = 1      {category: Initials; units: g/L}              [Initial glucose concentration]
  Gln = 1      {category: Initials; units: g/L}              [Initial glutamine concentration]
  Lac = 0.5    {category: Initials; units: g/L}              [Initial lactate concentration]
  Amm = 1.1    {category: Initials; units: g/L}              [Initial ammonium concentration]
  MAb = 0      {category: Initials; units: g/L}              [Initial antibody concentration]
  ATP = 0      {category: Initials; units: g/L}              [Initial ATP concentration]
  Vol = 2      {category: Initials; units: L; format: #0.00} [Initial volume] 

#parameters:
  mu_1_max        = 0.048     {caption:  mu_1_max ; category:  Parameters ; units:  1/h ; min: 0.00000001; max: 1500; step: 0.00001} [Maximum growth rate (exponential)]
  mu_2_max        = 0.012     {caption:  mu_2_max ; category:  Parameters ; units:  1/h ; min: 0.00000001; max: 1500; step: 0.00001} [Maximum growth rate (stationary)]
  K_Glc           = 1         {caption:  K_Glc ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 10; step: 0.0000001} [Monod constant of glucose]
  K_Gln           = 0.22      {caption:  K_Gln ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 10; step: 0.0000001} [Monod constant of glutamine]
  K_Lac           = 0.2       {caption:  K_Lac ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 0.3; step: 0.0000001} [Monod constant of lactate]
  K_Lac_i         = 150       {caption:  K_Lac_i ; category:  Parameters ; units: g/L ; min: 0.00000001; max: 30000.0; step: 1.0} [Constant of lactate inhibition]
  K_Amm           = 40        {category:  Parameters}
  k_mu            = 0.01      {category:  Parameters}
  k_d_max         = 0.01      {category:  Parameters}
  m_Glc           = 8e-6      {category:  Parameters; format: #0.0E00}
  m_Gln           = 0         {category:  Parameters}
  Y_X_Glc         = 15        {category:  Parameters}
  Y_X_Gln         = 190       {category:  Parameters}
  Y_Lac_Glc       = 4.0       {category:  Parameters}
  Y_X_Lac         = 7.500075  {category:  Parameters; min: 0.0001; max: 30}
  Y_Amm_Gln       = 1.1       {category:  Parameters}
  K_lysis         = 0.000001  {category:  Parameters}
  Glc_Cell_Th     = 0.2       {category:  Parameters; min: 0.001; max: 0.4}
  sigm_factor     = -10       {category:  Parameters; min: -100; max: -10}
  W_dry           = 2.5e-12   {category:  Constants; units: g/cell; format: #0.0E00}

#constants:
  W_dry_scale = 1e7

#tolerance: 0.0000001`;

export const feeding = new Map([
  [24, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [48, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [72, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [96, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [120, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [144, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [168, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [192, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
  [216, new Float32Array([0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])],
]);

/** Names of "main" feeding functions */
type FeedNames = {
    vol: string,
    glc: string,
    gln: string,
  };

  /** Default argument names*/
  enum ARG {
    START = '_t0',
    FINISH = '_t1',
    STEP = '_h',
  };

import * as DSL from '@datagrok/diff-grok';

const START_IDX = 0;
const FINISH_IDX = 1;
const ARG_COL_IDX = 0;

/** Pipeline creator for Dynamic Digital Twin */
export class DdtPipelineCreator extends DSL.PipelineCreator {
  private feeding: Map<number, Float32Array>;
  private ouputCode: string | null;
  private inputNames: string[];
  private feedNames: FeedNames;
  private funcNames: string[];
  private correctionFactor: number;
  private defaultCorrection: number;

  constructor(ivp: DSL.IVP, feeding: Map<number, Float32Array>, feedNames: FeedNames, funcNames: string[],
    correctionFactor: number, defaultCorrection: number) {
    super(ivp);
    this.feeding = feeding;
    this.ouputCode = DSL.getOutputCode(this.ivp);
    this.inputNames = ([ARG.START, ARG.FINISH, ARG.START] as string[]).concat(ivp.deqs.solutionNames);
    this.feedNames = feedNames;
    this.funcNames = funcNames;
    this.correctionFactor = correctionFactor;
    this.defaultCorrection = defaultCorrection;
  } // constructor

  /** Return simulation pipeline */
  public getPipeline(inputs: Float64Array): DSL.Pipeline {
    const VOL_IDX = this.inputNames.indexOf(this.feedNames.vol);
    const GLC_IDX = this.inputNames.indexOf(this.feedNames.glc);
    const GLN_IDX = this.inputNames.indexOf(this.feedNames.gln);

    const globalFinish = inputs[FINISH_IDX];

    //console.log(this.feeding);

    const feedingTimes: number[] = [];
    this.feeding.forEach((_, time) => {
      if (time <= globalFinish)
        feedingTimes.push(time);
    });
    feedingTimes.sort((a, b) => a - b);

    //console.log(feedingTimes);
    //console.log(globalFinish);

    const timesCount = feedingTimes.length;

    if (timesCount < 1) {
      return {
        wrappers: [
          {
            preproc: null,
            out: this.ouputCode,
            postproc: null,
          },
        ],
        out: null,
      };
    }

    let V_delta: number;
    let rawVals: Float32Array;

    const wrappers: DSL.Wrapper[] = [];
    let finish = 0;
    let preproc = '';
    let postproc = '';
    let inIdx = 0;
    let outIdx = 0;

    let lines: string[];

    for (let idx = 0; idx < timesCount; ++idx) {
      finish = feedingTimes[idx];

      preproc = `arguments[1][${FINISH_IDX}] = ${finish}; return arguments[1];`;

      rawVals = this.feeding.get(finish)!;

      V_delta = rawVals[VOL_IDX];

      lines = [
        `const V = arguments[1][${VOL_IDX}];`,
        'const lastIndex = arguments[0][0].length - 1;',
      ];

      //  ['X_T', 'X_V', 'Glc', 'Gln', 'Lac', 'Amm', 'MAb', 'ATP'].forEach((name) => {
      //    inputs[nameToIndexMap.get(name) ?? 0] = solution!.get(name, lastRowIdx) /(1 + V_delta/V);
      //  });
      this.funcNames.forEach((name) => {
        inIdx = this.inputNames.indexOf(name);
        if (inIdx < 0)
          throw new Error(`Inconsistent model: no '${name}' in inputs. Fix the model`);

        outIdx = this.outputNames.indexOf(name);
        if (outIdx < 0)
          throw new Error(`Inconsistent model: no '${name}' in outputs. Fix the model`);

        lines.push(`arguments[1][${inIdx}] = arguments[0][${outIdx}][lastIndex] / (1 + ${V_delta} / V);`);
      });

      // Update Glc
      // inputs[GLC_IDX] += rawVals[GLC_IDX] * V_delta / (V + V_delta);
      lines.push(
        `arguments[1][${GLC_IDX}] = arguments[1][${GLC_IDX}] + ${rawVals[GLC_IDX]} * ${V_delta} / (V + ${V_delta});`,
      );

      // Update Gln
      // inputs[GLN_IDX] += rawVals[GLN_IDX] * V_delta / (V + V_delta);
      lines.push(
        `arguments[1][${GLN_IDX}] = arguments[1][${GLN_IDX}] + ${rawVals[GLN_IDX]} * ${V_delta} / (V + ${V_delta});`,
      );

      // Update volume
      // V += V_delta;
      // inputs[VOL_IDX] = V;
      lines.push(`arguments[1][${VOL_IDX}] = arguments[1][${VOL_IDX}] + ${V_delta};`);

      //console.log('Feeding');
      //console.log(lines);

      postproc = `arguments[1][${START_IDX}] = ${finish};
        ${lines.join('\n')}
        return arguments[1];`;

      wrappers.push({
        preproc: preproc,
        out: this.ouputCode,
        postproc: postproc,
      });
    }

    //console.log(wrappers);

    if (globalFinish > finish) {
      wrappers.push({
        preproc: `arguments[1][${FINISH_IDX}] = ${globalFinish};return arguments[1];`,
        out: this.ouputCode,
        postproc: null,
      });
    }

    return {
      wrappers: wrappers,
      out: this.getArgumentCorrectionCode(),
    };
  } // getPipeline

  /** Return code for correction of an argument */
  private getArgumentCorrectionCode(): string {
    return `const arg = arguments[0][${ARG_COL_IDX}];
      const length = arg.length;
      const step = arg[1] - arg[0];
      const correction = Math.min(step * ${this.correctionFactor}, ${this.defaultCorrection});
      for (let i = 0; i < length - 1; ++i) {
        if (arg[i + 1] - arg[i] < step)
          arg[i] -= correction;
      }
      return arguments[0];`;
  } // getArgumentCorrectionCode
} // DdtPipelineCreator
