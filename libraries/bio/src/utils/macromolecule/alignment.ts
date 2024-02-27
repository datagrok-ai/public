import {GapSeqMonomer, ISeqMonomer, ISeqSplitted} from './types';
import {GapSymbols} from '../units-handler';
import {NOTATION} from './consts';

export type AlignmentOptions = {
  gapPenalty: number,
  matchScore: number,
  gapSymbol: ISeqMonomer,
  localAlignment: boolean,
}

export const defaultAlignmentOptions: AlignmentOptions = {
  gapPenalty: 1,
  matchScore: 1,
  gapSymbol: GapSymbols[NOTATION.FASTA],
  localAlignment: false,
};

export function alignSequencePair(template: ISeqSplitted, seq: ISeqSplitted,
  alignmentOptions: Partial<AlignmentOptions> = {},
): { seq1: string, seq2: string, seq1Splitted: ISeqMonomer[], seq2Splitted: ISeqMonomer[] } {
  const options = {...defaultAlignmentOptions, ...alignmentOptions};

  function compare(c1: ISeqMonomer, c2: ISeqMonomer) {
    return c1.canonical === c2.canonical ? options.matchScore : -options.matchScore;
  }

  const templateLen = template.length;
  const seqLen = seq.length;

  const matrix = Array(templateLen + 1).fill(0).map(() => Array(seqLen + 1).fill(0));

  let maxScore = -9999;
  let index = [1, 1];
  for (let i = 1; i < templateLen + 1; i++) {
    for (let j = 1; j < seqLen + 1; j++) {
      matrix[i][j] = Math.max(0,
        matrix[i - 1][j - 1] + compare(template[i - 1], seq[j - 1]),
        matrix[i - 1][j] - options.gapPenalty,
        matrix[i][j - 1] - options.gapPenalty);
      if (matrix[i][j] >= maxScore) {
        maxScore = matrix[i][j];
        index = [i, j];
      }
    }
  }

  let i = options.localAlignment ? index[0] : templateLen;
  let j = options.localAlignment ? index[1] : seqLen;
  let alignedTemplate = new Array<ISeqMonomer>(0);
  let alignedSeq = new Array<ISeqMonomer>(0);

  while (i > 0 && j > 0) {
    const mx = Math.max(matrix[i - 1][j - 1], matrix[i - 1][j], matrix[i][j - 1]);

    if (matrix[i][j] == matrix[i - 1][j - 1] + options.matchScore && mx == matrix[i - 1][j - 1]) {
      alignedTemplate.push(template[i - 1]);
      alignedSeq.push(seq[j - 1]);
      i -= 1;
      j -= 1;
    } else {
      if (matrix[i][j] == matrix[i - 1][j] - options.gapPenalty) {
        alignedSeq.push(options.gapSymbol);
        alignedTemplate.push(template[i - 1]);
        i -= 1;
      } else {
        if (matrix[i][j] == matrix[i][j - 1] - options.gapPenalty) {
          alignedTemplate.push(options.gapSymbol);
          alignedSeq.push(seq[j - 1]);
          j -= 1;
        } else {
          alignedTemplate.push(template[i - 1]);
          alignedSeq.push(seq[j - 1]);
          i -= 1;
          j -= 1;
        }
      }
    }
  }
  alignedTemplate = [...Array.from(template).splice(0, i), ...alignedTemplate.reverse(),
    ...(options.localAlignment ? Array.from(template).splice(index[0], templateLen) : [])];
  alignedSeq = [...Array.from(seq).splice(0, j), ...alignedSeq.reverse(),
    ...(options.localAlignment ? Array.from(seq).splice(index[1], seqLen) : [])];

  const templateStart = i;
  const seqStart = j;

  if (templateStart > seqStart)
    alignedSeq = [...new Array(templateStart - seqStart).fill(options.gapSymbol.valueOf()), ...alignedSeq];
  else
    alignedTemplate = [...new Array(seqStart - templateStart).fill(options.gapSymbol.valueOf()), ...alignedTemplate];


  if (alignedSeq.length > alignedTemplate.length)
    alignedTemplate.push(...new Array(alignedSeq.length - alignedTemplate.length).fill(options.gapSymbol.valueOf()));
  else
    alignedSeq.push(...new Array(alignedTemplate.length - alignedSeq.length).fill(options.gapSymbol.valueOf()));

  return {
    seq1: alignedTemplate.join(''), seq2: alignedSeq.join(''),
    seq1Splitted: alignedTemplate, seq2Splitted: alignedSeq
  };
}
