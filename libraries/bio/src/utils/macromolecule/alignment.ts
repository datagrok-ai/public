import wu from 'wu';

import {ISeqSplitted} from './types';
import {GapOriginals, NOTATION} from './consts';
import {StringListSeqSplitted} from './utils';

export type AlignmentOptions = {
  gapPenalty: number,
  matchScore: number,
  gapSymbol: string,
  localAlignment: boolean,
}

export const defaultAlignmentOptions: AlignmentOptions = {
  gapPenalty: 1,
  matchScore: 1,
  gapSymbol: GapOriginals[NOTATION.FASTA],
  localAlignment: false,
};

export function alignSequencePair(template: ISeqSplitted, seq: ISeqSplitted,
  alignmentOptions: Partial<AlignmentOptions> = {},
): { seq1: string, seq2: string, seq1Splitted: ISeqSplitted, seq2Splitted: ISeqSplitted } {
  const options = {...defaultAlignmentOptions, ...alignmentOptions};

  function compare(cm1: string, cm2: string) {
    return cm1 === cm2 ? options.matchScore : -options.matchScore;
  }

  const templateLen = template.length;
  const seqLen = seq.length;

  const matrix = Array(templateLen + 1).fill(0).map(() => Array(seqLen + 1).fill(0));

  let maxScore = -9999;
  let index = [1, 1];
  for (let i = 1; i < templateLen + 1; i++) {
    for (let j = 1; j < seqLen + 1; j++) {
      matrix[i][j] = Math.max(0,
        matrix[i - 1][j - 1] + compare(template.getCanonical(i - 1), seq.getCanonical(j - 1)),
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
  let alignedTemplate = new Array<string>(0);
  let alignedSeq = new Array<string>(0);

  while (i > 0 && j > 0) {
    const mx = Math.max(matrix[i - 1][j - 1], matrix[i - 1][j], matrix[i][j - 1]);

    if (matrix[i][j] == matrix[i - 1][j - 1] + options.matchScore && mx == matrix[i - 1][j - 1]) {
      alignedTemplate.push(template.getCanonical(i - 1));
      alignedSeq.push(seq.getCanonical(j - 1));
      i -= 1;
      j -= 1;
    } else {
      if (matrix[i][j] == matrix[i - 1][j] - options.gapPenalty) {
        alignedSeq.push(options.gapSymbol);
        alignedTemplate.push(template.getCanonical(i - 1));
        i -= 1;
      } else {
        if (matrix[i][j] == matrix[i][j - 1] - options.gapPenalty) {
          alignedTemplate.push(options.gapSymbol);
          alignedSeq.push(seq.getCanonical(j - 1));
          j -= 1;
        } else {
          alignedTemplate.push(template.getCanonical(i - 1));
          alignedSeq.push(seq.getCanonical(j - 1));
          i -= 1;
          j -= 1;
        }
      }
    }
  }
  alignedTemplate = [...wu.count(0).take(i).map((pos) => template.getCanonical(pos)), ...alignedTemplate.reverse(),
    ...(options.localAlignment ? wu.count(index[0]).take(templateLen).map((pos) => template.getCanonical(pos)) : [])];
  alignedSeq = [...wu.count(0).take(j).map((pos) => seq.getCanonical(pos)), ...alignedSeq.reverse(),
    ...(options.localAlignment ? wu.count(index[1]).take(seqLen).map((pos) => seq.getCanonical(pos)) : [])];

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
    seq1: alignedTemplate.join(''),
    seq2: alignedSeq.join(''),
    seq1Splitted: new StringListSeqSplitted(alignedTemplate, options.gapSymbol),
    seq2Splitted: new StringListSeqSplitted(alignedSeq, options.gapSymbol)
  };
}
