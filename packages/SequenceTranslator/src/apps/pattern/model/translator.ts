export function applyPatternToRawSequence(
  rawNucleotideSequence: string,
  modifications: string[],
  ptoFlags: boolean[]
): string {
  const rawNucleotides = rawNucleotideSequence.split('');
  const modifiedNucleotides = rawNucleotides.map((nucleotide, i) => {
    const modifiedNucleotide = getModifiedNucleotide(nucleotide, modifications[i]);
    return modifiedNucleotide;
  });

  const modificationsWithPTOLinkages = getModificationsWithPTOLinkages(modifiedNucleotides, ptoFlags);

  return modificationsWithPTOLinkages.join('');
}

function getModifiedNucleotide(nucleotide: string, modification: string): string {
  return `[${nucleotide}]`;
}

function getPhosphorothioateLinkageSymbol(): string {
  return 'ps';
}

function getModificationsWithPTOLinkages(
  modifiedNucleotides: string[],
  ptoFlags: boolean[]
): string[] {
  const modificationsWithPTOLinkages = new Array<string>(
    modifiedNucleotides.length + ptoFlags.filter((flag) => flag).length
  );

  const ptoLinkage = getPhosphorothioateLinkageSymbol();
  let idxShift = 0;

  if (ptoFlags[0]) {
    modificationsWithPTOLinkages[0] = ptoLinkage;
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

