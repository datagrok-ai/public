/**
 * Library of molecular weights.
 * @link https://worldwide.promega.com/resources/tools/amino-acid-chart-amino-acid-structure */
const _lib = [
  {
    'Name': 'Alanine',
    'Three-letter': 'Ala',
    'One-letter': 'A',
    'Weight': '89Da',
  },
  {
    'Name': 'Arginine',
    'Three-letter': 'Arg',
    'One-letter': 'R',
    'Weight': '174Da',
  },
  {
    'Name': 'Asparagine',
    'Three-letter': 'Asn',
    'One-letter': 'N',
    'Weight': '132Da',
  },
  {
    'Name': 'Aspartic acid',
    'Three-letter': 'Asp',
    'One-letter': 'D',
    'Weight': '133Da',
  },
  {
    'Name': 'Asparagine or aspartic acid',
    'Three-letter': 'Asx',
    'One-letter': 'B',
    'Weight': '133Da',
  },
  {
    'Name': 'Cysteine',
    'Three-letter': 'Cys',
    'One-letter': 'C',
    'Weight': '121Da',
  },
  {
    'Name': 'Glutamine',
    'Three-letter': 'Gln',
    'One-letter': 'Q',
    'Weight': '146Da',
  },
  {
    'Name': 'Glutamic acid',
    'Three-letter': 'Glu',
    'One-letter': 'E',
    'Weight': '147Da',
  },
  {
    'Name': 'Glutamine or glutamic acid',
    'Three-letter': 'Glx',
    'One-letter': 'Z',
    'Weight': '147Da',
  },
  {
    'Name': 'Glycine',
    'Three-letter': 'Gly',
    'One-letter': 'G',
    'Weight': '75Da',
  },
  {
    'Name': 'Histidine',
    'Three-letter': 'His',
    'One-letter': 'H',
    'Weight': '155Da',
  },
  {
    'Name': 'Isoleucine',
    'Three-letter': 'Ile',
    'One-letter': 'I',
    'Weight': '131Da',
  },
  {
    'Name': 'Leucine',
    'Three-letter': 'Leu',
    'One-letter': 'L',
    'Weight': '131Da',
  },
  {
    'Name': 'Lysine',
    'Three-letter': 'Lys',
    'One-letter': 'K',
    'Weight': '146Da',
  },
  {
    'Name': 'Methionine',
    'Three-letter': 'Met',
    'One-letter': 'M',
    'Weight': '149Da',
  },
  {
    'Name': 'Phenylalanine',
    'Three-letter': 'Phe',
    'One-letter': 'F',
    'Weight': '165Da',
  },
  {
    'Name': 'Proline',
    'Three-letter': 'Pro',
    'One-letter': 'P',
    'Weight': '115Da',
  },
  {
    'Name': 'Serine',
    'Three-letter': 'Ser',
    'One-letter': 'S',
    'Weight': '105Da',
  },
  {
    'Name': 'Threonine',
    'Three-letter': 'Thr',
    'One-letter': 'T',
    'Weight': '119Da',
  },
  {
    'Name': 'Tryptophan',
    'Three-letter': 'Trp',
    'One-letter': 'W',
    'Weight': '204Da',
  },
  {
    'Name': 'Tyrosine',
    'Three-letter': 'Tyr',
    'One-letter': 'Y',
    'Weight': '181Da',
  },
  {
    'Name': 'Valine',
    'Three-letter': 'Val',
    'One-letter': 'V',
    'Weight': '117Da',
  },
];

/**
 * Selection of the molecular weights of amino acid residues.
 * @type {*} */
const weightsLib : {[name: string]: number} = {};

// Create a dictionary linking one-letter code with the corresponding residues weight.
for (const d of _lib)
  weightsLib[d['One-letter']] = parseFloat(d.Weight.substring(0, d.Weight.length-2));


/**
 * Calculates molecular weight of the given peptide in daltons.
 *
 * @export
 * @param {string} sequence Peptide sequence.
 * @return {number} Molecular weight in Da.
 */
export function getSequenceMolecularWeight(sequence: string): number {
  let sum = 0;

  if (sequence.startsWith('NH2')) {
    sum += 16.02;
    sequence = sequence.substring(3);
  }

  if (sequence.endsWith('COOH')) {
    sum += 45.02;
    sequence = sequence.substring(0, sequence.length-4);
  }

  for (const i of sequence) {
    if (i in weightsLib)
      sum += weightsLib[i];
  }
  return sum;
}
