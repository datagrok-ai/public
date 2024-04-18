import {TERMINI, TERMINUS} from './const';
import {PATTERN_APP_DATA} from '../../common/model/data-loader/json-loader';

export function applyPatternToRawSequence(
  rawNucleotideSequence: string,
  modifications: string[],
  ptoFlags: boolean[],
  terminalModifications: Record<TERMINUS, string>
): string {
  const rawNucleotides = rawNucleotideSequence.split('');

  const modifiedNucleotides = rawNucleotides.map((nucleotide, i) => {
    const modifiedNucleotide = getModifiedNucleotide(nucleotide, modifications[i]);
    return modifiedNucleotide;
  });

  const modificationsWithPTOLinkages = getModificationsWithPTOLinkages(
    modifiedNucleotides, ptoFlags, terminalModifications
  );

  return modificationsWithPTOLinkages.join('');
}

function getModifiedNucleotide(nucleotide: string, modification: string): string {
  const substitution = PATTERN_APP_DATA[modification]['substitution'];
  return nucleotide.replace(/([AGCTU])/, substitution);
}

function getPhosphorothioateLinkageSymbol(): string {
  return 'ps';
}

function getModificationsWithPTOLinkages(
  modifiedNucleotides: string[],
  ptoFlags: boolean[],
  terminalModifications: Record<TERMINUS, string>
): string[] {
  const modificationsWithPTOLinkages = new Array<string>(
    modifiedNucleotides.length + ptoFlags.filter((flag) => flag).length + TERMINI.length
  );

  const ptoLinkage = getPhosphorothioateLinkageSymbol();

  modificationsWithPTOLinkages[0] = terminalModifications[TERMINUS.FIVE_PRIME];
  modificationsWithPTOLinkages[modificationsWithPTOLinkages.length - 1] = terminalModifications[TERMINUS.THREE_PRIME];

  let idxShift = 1;

  if (ptoFlags[0]) {
    modificationsWithPTOLinkages[idxShift] = ptoLinkage;
    idxShift++;
  }

  modifiedNucleotides.forEach((nucleotide, i) => {
    modificationsWithPTOLinkages[i + idxShift] = nucleotide;
    if (ptoFlags[i + 1]) {
      modificationsWithPTOLinkages[i + idxShift + 1] = ptoLinkage;
      idxShift++;
    }
  });

  return modificationsWithPTOLinkages;
}
