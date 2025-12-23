import { Sample, DescriptorStats } from './types';
import {
  mean,
  std,
  gaussian,
  welchTTest,
  pearson,
} from './stats';

export interface PMPOOptions {
  pCutoff?: number;
  corrCutoff?: number;
  sigmoidalCorrection?: boolean;
}

export class pMPOBuilder {
  descriptors: DescriptorStats[] = [];

  constructor(
    data: Sample[],
    public labelKey: string,
    public goodLabel: any,
    options: PMPOOptions = {}
  ) {
    this.pCutoff = options.pCutoff ?? 0.01;
    this.corrCutoff = options.corrCutoff ?? 0.75;
    this.sigmoidalCorrection =
      options.sigmoidalCorrection ?? true;

    this.buildModel(data);
  }

  pCutoff: number;
  corrCutoff: number;
  sigmoidalCorrection: boolean;

  buildModel(data: Sample[]) {
    const keys = Object.keys(data[0]).filter(
      (k) => k !== this.labelKey
    );

    /* ---------- Step 1: t-test filtering ---------- */

    for (const key of keys) {
      const good = data
        .filter((d) => d[this.labelKey] === this.goodLabel)
        .map((d) => d[key])
        .filter((v): v is number => typeof v === 'number');

      const bad = data
        .filter((d) => d[this.labelKey] !== this.goodLabel)
        .map((d) => d[key])
        .filter((v): v is number => typeof v === 'number');

      if (good.length < 3 || bad.length < 3) continue;

      const mg = mean(good);
      const mb = mean(bad);
      const sg = std(good, mg);
      const sb = std(bad, mb);

      const { t, p } = welchTTest(good, bad);
      if (p > this.pCutoff) continue;

      this.descriptors.push({
        name: key,
        meanGood: mg,
        meanBad: mb,
        stdGood: sg,
        stdBad: sb,
        weight: Math.abs(mg - mb),
        t,
        p,
      });
    }

    /* ---------- Step 2: correlation filtering ---------- */

    this.filterCorrelated(data);
  }

  private filterCorrelated(data: Sample[]) {
    const kept: DescriptorStats[] = [];

    for (const d of this.descriptors.sort(
      (a, b) => a.p - b.p
    )) {
      let correlated = false;

      for (const k of kept) {
        const x = data
          .map((r) => r[d.name])
          .filter((v): v is number => typeof v === 'number');

        const y = data
          .map((r) => r[k.name])
          .filter((v): v is number => typeof v === 'number');

        if (x.length !== y.length) continue;

        const rxy = Math.abs(pearson(x, y));
        if (rxy >= this.corrCutoff) {
          correlated = true;
          break;
        }
      }

      if (!correlated) kept.push(d);
    }

    this.descriptors = kept;
  }

  /* ---------- Scoring ---------- */

  score(sample: Sample): number {
    let s = 0;

    for (const d of this.descriptors) {
      const x = sample[d.name] as number;
      if (typeof x !== 'number') continue;

      const g = gaussian(x, d.meanGood, d.stdGood || 1);
      const b = gaussian(x, d.meanBad, d.stdBad || 1);

      let term = d.weight * (g / (g + b));

      if (this.sigmoidalCorrection) {
        term *= 1 / (1 + Math.exp(-(x - d.meanGood)));
      }

      s += term;
    }

    return s;
  }

  get statistics(): DescriptorStats[] {
    return [...this.descriptors];
  }
}
