import weighted from 'weighted';

export class SequenceGenerator {
  static maxPopulationLength = 4;
  alphabet = 'ACDEFGHIKLMNPQRSTVWY-';
  length: number;
  alphas: string[];

  constructor(length) {
    this.length = length;
  }

  genList(): string[] {
    const cutoff = Math.random();
    const alphabet = Array.from(this.alphabet);
    let alpha = alphabet.filter((_) => Math.random() > cutoff);

    if (alpha.length == 0)
      alpha = [weighted.select(alphabet)];

    const weights = alpha.map((_) => Math.random());
    const sum = weights.reduce((prev, curr) => prev + curr);
    const portions = weights.map((v) => v / sum);
    const seqs: string[] = [];

    for (let i = 0; i < this.length; ++i)
      seqs.push(weighted.select(alpha, portions));

    return seqs;
  }
}
