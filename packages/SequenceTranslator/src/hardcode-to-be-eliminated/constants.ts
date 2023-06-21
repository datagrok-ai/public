const RNA_COLOR = 'rgb(255,230,153)';
const INVABASIC_COLOR = 'rgb(203,119,211)';
export const AXOLABS_MAP:
  {[index: string]: {fullName: string, symbols: [string, string, string, string], color: string}} =
{
  'RNA': {
    fullName: 'RNA nucleotides',
    symbols: ['A', 'C', 'G', 'U'],
    color: RNA_COLOR,
  },
  'DNA': {
    fullName: 'DNA nucleotides',
    symbols: ['dA', 'dC', 'dG', 'dT'],
    color: 'rgb(197,224,180)',
  },
  "2\"-Fluoro": {
    fullName: '2\'-Fluoro nucleotides',
    symbols: ['Af', 'Cf', 'Gf', 'Uf'],
    color: 'rgb(68,114,196)',
  },
  "2\"-O-Methyl": {
    fullName: '2\'-O-Methyl nucleotides',
    symbols: ['a', 'c', 'g', 'u'],
    color: 'rgb(166,166,166)',
  },
  "2\"-O-MOE": {
    fullName: '2\'-O-MOE nucleotides (including 5-Methyl C)',
    symbols: ['Am', 'Cm', 'Gm', 'Tm'],
    color: 'rgb(112,48,160)',
  },
  'GNA': {
    fullName: 'Glycol nucleic acid',
    symbols: ['(GNA-A)', '(GNA-C)', '(GNA-G)', '(GNA-T)'],
    color: 'rgb(255,192,0)',
  },
  'LNA': {
    fullName: 'Locked nucleic acid (including 5-Methyl C)',
    symbols: ['Ab', 'Cb', 'Gb', 'Tb'],
    color: 'rgb(54,229,238)',
  },
  'UNA': {
    fullName: 'Unlocked nucleotides',
    symbols: ['Ao', 'Co', 'Go', 'Uo'],
    color: 'rgb(255,192,0)',
  },
  'A': {
    fullName: 'Adenosine',
    symbols: ['a', 'a', 'a', 'a'],
    color: RNA_COLOR,
  },
  'C': {
    fullName: 'Cytidine',
    symbols: ['c', 'c', 'c', 'c'],
    color: RNA_COLOR,
  },
  'G': {
    fullName: 'Guanosine',
    symbols: ['g', 'g', 'g', 'g'],
    color: RNA_COLOR,
  },
  'U': {
    fullName: 'Uridine',
    symbols: ['u', 'u', 'u', 'u'],
    color: RNA_COLOR,
  },
  'X-New': {
    fullName: '',
    symbols: ['X', 'X', 'X', 'X'],
    color: 'rgb(108,0,0)',
  },
  'Y-New': {
    fullName: '',
    symbols: ['Y', 'Y', 'Y', 'Y'],
    color: 'rgb(210,146,146)',
  },
  'Z-New': {
    fullName: '',
    symbols: ['Z', 'Z', 'Z', 'Z'],
    color: 'rgb(155,108,132)',
  },
  'InvAbasic': {
    fullName: 'Inverted abasic capped',
    symbols: ['(invabasic)', '(invabasic)', '(invabasic)', '(invabasic)'],
    color: INVABASIC_COLOR,
  },
  "5\"-vinylps": {
    fullName: '5\'-vinylphosphonate-2\'-OMe-uridine',
    symbols: ['(vinu)', '(vinu)', '(vinu)', '(vinu)'],
    color: 'rgb(0,0,139)',
  },
  'InvAbasic(o)': {
    fullName: 'Inverted abasic capped (overhang)',
    symbols: ['(invabasic)', '(invabasic)', '(invabasic)', '(invabasic)'],
    color: INVABASIC_COLOR,
  },
  "2\"-OMe-U(o)": {
    fullName: 'Nucleotide Uridine with 2\'O-Methyl protection (overhang)',
    symbols: ['mU', 'mU', 'mU', 'mU'],
    color: 'rgb(65,233,80)',
  },
};
