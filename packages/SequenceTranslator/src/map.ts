export const SYNTHESIZERS = {
  RAW_NUCLEOTIDES: 'Raw Nucleotides',
  BIOSPRING: 'BioSpring Codes',
  GCRS: 'Janssen GCRS Codes',
  AXOLABS: 'Axolabs Codes',
  MERMADE_12: 'Mermade 12',
};
export const TECHNOLOGIES = {
  DNA: 'DNA',
  RNA: 'RNA',
  ASO_GAPMERS: 'For ASO Gapmers',
  SI_RNA: 'For 2\'-OMe and 2\'-F modified siRNA',
};
// interface CODES {
// }
export const MODIFICATIONS: {[index: string]: {left: string, right: string}} = {
  '(invabasic)': {
    left: 'O[C@@H]1C[C@@H]O[C@H]1CO',
    right: 'O[C@@H]1C[C@@H]O[C@H]1CO',
  },
  '(GalNAc-2-JNJ)': {
    left: 'C(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)NC(=O)CCCC(=O)NCC(O)CO',
    right: 'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)'+
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)',
  },
};
export const stadardPhosphateLinkSmiles = 'OP(=O)(O)O';
export const map: {[synthesizer: string]:
  {[technology: string]: {[code: string]:
    {'name': string, 'weight': number, 'normalized': string, 'SMILES': string}}}} = {
      'Raw Nucleotides': {
        'DNA': {
          'A': {
            'name': 'Adenine',
            'weight': 313.21,
            'normalized': 'dA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O',
          },
          'T': {
            'name': 'Tyrosine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
          'G': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O',
          },
          'C': {
            'name': 'Cytosine',
            'weight': 289.18,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O',
          },
        },
        'RNA': {
          'A': {
            'name': 'Adenine',
            'weight': 313.21,
            'normalized': 'dA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O',
          },
          'U': {
            'name': 'Uracil',
            'weight': 306.17,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](O)[C@@H]1O',
          },
          'G': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O',
          },
          'C': {
            'name': 'Cytosine',
            'weight': 289.18,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O',
          },
        },
      },
      'BioSpring Codes': {
        'For ASO Gapmers': {
          '5': {
            'name': '2\'MOE-5Me-rU',
            'weight': 378.27,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))[C@H](OCCOC)[C@@H]1O',
          },
          '6': {
            'name': '2\'MOE-rA',
            'weight': 387.29,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
          },
          '7': {
            'name': '2\'MOE-5Me-rC',
            'weight': 377.29,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O',
          },
          '8': {
            'name': '2\'MOE-rG',
            'weight': 403.28,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
          },
          '9': {
            'name': '5-Methyl-dC',
            'weight': 303.21,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O',
          },
          '*': {
            'name': 'ps linkage',
            'weight': 16.07,
            'normalized': '',
            'SMILES': 'OP(=O)(S)O',
          },
          'A': {
            'name': 'Adenine',
            'weight': 313.21,
            'normalized': 'dA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O',
          },
          'C': {
            'name': 'Cytosine',
            'weight': 289.18,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O',
          },
          'G': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O',
          },
          'T': {
            'name': 'Tyrosine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
        },
        'For 2\'-OMe and 2\'-F modified siRNA': {
          '1': {
            'name': '2\'-fluoro-U',
            'weight': 308.16,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',
          },
          '2': {
            'name': '2\'-fluoro-A',
            'weight': 331.2,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          '3': {
            'name': '2\'-fluoro-C',
            'weight': 307.18,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',
          },
          '4': {
            'name': '2\'-fluoro-G',
            'weight': 347.19,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          '5': {
            'name': '2\'OMe-rU',
            'weight': 320.2,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',
          },
          '6': {
            'name': '2\'OMe-rA',
            'weight': 343.24,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          '7': {
            'name': '2\'OMe-rC',
            'weight': 319.21,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',
          },
          '8': {
            'name': '2\'OMe-rG',
            'weight': 359.24,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          '*': {
            'name': 'ps linkage',
            'weight': 16.07,
            'normalized': '',
            'SMILES': 'OP(=O)(S)O',
          },
        },
      },
      'Axolabs Codes': {
        'For 2\'-OMe and 2\'-F modified siRNA': {
          'Uf': {
            'name': '2\'-fluoro-U',
            'weight': 308.16,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',
          },
          'Af': {
            'name': '2\'-fluoro-A',
            'weight': 331.2,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'Cf': {
            'name': '2\'-fluoro-C',
            'weight': 307.18,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',
          },
          'Gf': {
            'name': '2\'-fluoro-G',
            'weight': 347.19,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'u': {
            'name': '2\'OMe-rU',
            'weight': 320.2,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'a': {
            'name': '2\'OMe-rA',
            'weight': 343.24,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          'c': {
            'name': '2\'OMe-rC',
            'weight': 319.21,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'g': {
            'name': '2\'OMe-rG',
            'weight': 359.,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          's': {
            'name': 'ps linkage',
            'weight': 16.07,
            'normalized': '',
            'SMILES': 'OP(=O)(S)O',
          },
        },
      },
      'Janssen GCRS Codes': {
        'For ASO Gapmers': {
          'moeT': {
            'name': '2\'MOE-5Me-rU',
            'weight': 378.27,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))[C@H](OCCOC)[C@@H]1O',
          },
          'moeA': {
            'name': '2\'MOE-rA',
            'weight': 387.29,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
          },
          'moe5mC': {
            'name': '2\'MOE-5Me-rC',
            'weight': 377.29,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O',
          },
          '(5m)moeC': {
            'name': '2\'MOE-5Me-rC',
            'weight': 377.29,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O',
          },
          'moeG': {
            'name': '2\'MOE-rG',
            'weight': 403.28,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
          },
          '5mC': {
            'name': '5-Methyl-dC',
            'weight': 303.28,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O',
          },
          '(5m)C': {
            'name': '5-Methyl-dC',
            'weight': 303.28,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O',
          },
          'ps': {
            'name': 'ps linkage',
            'weight': 16.07,
            'normalized': '',
            'SMILES': 'OP(=O)(S)O',
          },
          'A': {
            'name': 'Adenine',
            'weight': 313.21,
            'normalized': 'dA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O',
          },
          'dA': {
            'name': 'Adenine',
            'weight': 313.21,
            'normalized': 'dA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O',
          },
          'C': {
            'name': 'Cytosine',
            'weight': 289.18,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O',
          },
          'dC': {
            'name': 'Cytosine',
            'weight': 289.18,
            'normalized': 'dC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O',
          },
          'G': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O',
          },
          'dG': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O',
          },
          'T': {
            'name': 'Tyrosine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
          'dT': {
            'name': 'Tyrosine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
          'rA': {
            'name': 'Adenine',
            'weight': 329.21,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](O)[C@@H]1O',
          },
          'rC': {
            'name': 'Cytosine',
            'weight': 305.18,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](O)[C@@H]1O',
          },
          'rG': {
            'name': 'Guanine',
            'weight': 345.21,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](O)[C@@H]1O',
          },
          'rU': {
            'name': 'Uracil',
            'weight': 306.17,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](O)[C@@H]1O',
          },
        },
        'For 2\'-OMe and 2\'-F modified siRNA': {
          'fU': {
            'name': '2\'-fluoro-U',
            'weight': 308.16,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',
          },
          'fA': {
            'name': '2\'-fluoro-A',
            'weight': 331.2,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'fC': {
            'name': '2\'-fluoro-C',
            'weight': 307.18,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',
          },
          'fG': {
            'name': '2\'-fluoro-G',
            'weight': 347.19,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'mU': {
            'name': '2\'OMe-rU',
            'weight': 320.2,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'mA': {
            'name': '2\'OMe-rA',
            'weight': 343.24,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          'mC': {
            'name': '2\'OMe-rC',
            'weight': 319.21,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'mG': {
            'name': '2\'OMe-rG',
            'weight': 359.24,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
        },
      },
      'Mermade 12': {
        'For 2\'-OMe and 2\'-F modified siRNA': {
          'e': {
            'name': '2\'OMe-rA-ps',
            'weight': 359.31,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)O',
          },
          'h': {
            'name': '2\'OMe-rU-ps',
            'weight': 336.27,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)O',
          },
          'g': {
            'name': '2\'OMe-rG-ps',
            'weight': 375.31,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)O',
          },
          'f': {
            'name': '2\'OMe-rC-ps',
            'weight': 335.28,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)O',
          },
          'i': {
            'name': '2\'-fluoro-A-ps',
            'weight': 347.27,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)O',
          },
          'l': {
            'name': '2\'-fluoro-U-ps',
            'weight': 324.23,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)O',
          },
          'k': {
            'name': '2\'-fluoro-G-ps',
            'weight': 363.26,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)O',
          },
          'j': {
            'name': '2\'-fluoro-C-ps',
            'weight': 323.25,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)O',
          },
          'L': {
            'name': '2\'-fluoro-U',
            'weight': 308.16,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',
          },
          'I': {
            'name': '2\'-fluoro-A',
            'weight': 331.2,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'J': {
            'name': '2\'-fluoro-C',
            'weight': 307.18,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',
          },
          'K': {
            'name': '2\'-fluoro-G',
            'weight': 347.19,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O',
          },
          'H': {
            'name': '2\'OMe-rU',
            'weight': 320.2,
            'normalized': 'rU',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'E': {
            'name': '2\'OMe-rA',
            'weight': 343.24,
            'normalized': 'rA',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
          'F': {
            'name': '2\'OMe-rC',
            'weight': 319.21,
            'normalized': 'rC',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',
          },
          'G': {
            'name': '2\'OMe-rG',
            'weight': 359.24,
            'normalized': 'rG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',
          },
        },
      },
    };
