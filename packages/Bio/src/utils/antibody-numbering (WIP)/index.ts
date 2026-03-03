/**
 * antibody-numbering - Browser-compatible antibody sequence numbering library.
 *
 * Assigns antibody numbering (IMGT/Kabat/Chothia/AHo) to amino acid sequences
 * using profile-based alignment against consensus sequences.
 *
 * @example
 * ```ts
 * import { numberSequence } from 'antibody-numbering';
 *
 * const result = numberSequence(
 *   'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR',
 *   'imgt'
 * );
 *
 * console.log(result.chainType);        // 'Heavy'
 * console.log(result.numberingDetail);  // [{ position: '1', aa: 'Q' }, ...]
 * console.log(result.annotations);      // [{ name: 'FR1', ... }, ...]
 * ```
 */

export { numberSequence, numberSequences, extractSequence } from './annotator';
export type {
  Scheme,
  ChainType,
  ChainGroup,
  NumberingResult,
  NumberingEntry,
  RegionAnnotation,
  AlignmentResult,
} from './types';
