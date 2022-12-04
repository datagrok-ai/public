import * as DG from 'datagrok-api/dg';
import {getAllCodesOfSynthesizer} from './sequence-codes-tools';
import {differenceOfTwoArrays} from '../helpers';

export const DELIMITER = ';';
export const NUCLEOTIDES = ['A', 'G', 'C', 'U', 'T'];
export const SYNTHESIZERS = {
  RAW_NUCLEOTIDES: 'Raw Nucleotides',
  BIOSPRING: 'BioSpring Codes',
  GCRS: 'Janssen GCRS Codes',
  AXOLABS: 'Axolabs Codes',
  MERMADE_12: 'Mermade 12',
  LCMS: 'LCMS',
};
export const TECHNOLOGIES = {
  DNA: 'DNA',
  RNA: 'RNA',
  ASO_GAPMERS: 'For ASO Gapmers',
  SI_RNA: 'For 2\'-OMe and 2\'-F modified siRNA',
};
export const MODIFICATIONS: {[index: string]: {molecularWeight: number, left: string, right: string}} = {
  '(invabasic)': {
    molecularWeight: 118.13,
    left: 'O[C@@H]1C[C@@H]O[C@H]1CO',
    right: 'O[C@@H]1C[C@@H]O[C@H]1CO',
  },
  '(GalNAc-2-JNJ)': {
    molecularWeight: 1273.3,
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
            'name': 'Thymine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
          'G': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C[C@@H]1O',
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
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C[C@@H]1O',
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
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C[C@@H]1O',
          },
          'T': {
            'name': 'Thymine',
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
          's': {
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
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C[C@@H]1O',
          },
          'dG': {
            'name': 'Guanine',
            'weight': 329.21,
            'normalized': 'dG',
            'SMILES': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C[C@@H]1O',
          },
          'T': {
            'name': 'Thymine',
            'weight': 304.2,
            'normalized': 'dT',
            'SMILES': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O',
          },
          'dT': {
            'name': 'Thymine',
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
        'Others': {},
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

export const lcmsToGcrs = `LCMS, GCRS
A, A
C, C
/5mC/, (5m)C
G, G
T, T
rA, rA
rC, rC
rG, rG
rU, rU
mA, mA
mC, mC
/5mmC/, (5m)mC
mG, mG
mU, mU
fA, fA
fC, fC
/5mfC/, (5m)fC
fG, fG
fU, fU
/afA/, afA
/afC/, afC
/afG/, afG
/afU/, afU
+A, lna A
+C, lna C
+G, lna G
+T, lna T
/moeA/, moeA
/moeC/, moeC
/5mmoeC/, (5m)moeC
/moeG/, moeG
/moeT/, moeT
/moeU/, moeU
/xA/, Anp
/xC/, Cnp
/x5mC/, (5m)Cnp
/xG/, Gnp
/xT/, Tnp
/xrA/, rAnp
/xrC/, rCnp
/xrG/, rGnp
/xrU/, rUnp
/xmA/, mAnp
/xmC/, mCnp
/x5mmC/, (5m)mCnp
/xmG/, mGnp
/xmU/, mUnp
/xfA/, fAnp
/xfC/, fCnp
/xfG/, fGnp
/xfT/, fTnp
/xfU/, fUnp
/xafA/, afAnp
/xafC/, afCnp
/xafG/, afGnp
/xafU/, afUnp
/xeA/, eAnp
/xeC/, eCnp
/xeG/, eGnp
/xeU/, eUnp
/xmoeA/, moeAnp
/xmoeC/, moeCnp
/x5mmoeC/, (5m)moeCnp
/xmoeG/, moeGnp
/xmoeU/, moeUnp
/UNA-A/, (UNA-A)
/UNA-C/, (UNA-C)
/UNA-G/, (UNA-G)
/UNA-T/, (UNA-T)
/UNA-U/, (UNA-U)
/GNA-A/, (GNA-A)
/GNA-C/, (GNA-C)
/GNA-G/, (GNA-G)
/GNA-T/, (GNA-T)
/GNA-U/, (GNA-U)
/5CholTEG/, (5-CholTEG)
/3CholTEG/, (TEGChol-3)
/Toco/, Toco
/Palm/, Palm
/GalNAc/, GalNAc
/GalNAc2/, GalNAc2
/GalNAc3/, GalNAc3
/GalNAc6/, GalNAc6
/GalNAc7/, GalNAc7
/GalNAc9/, GalNAc9
/GalNAc14/, GalNAc14
/NAG37/, NAG37
/HEG/, (HEG)
/TEG/, (TEG)
/AmmC6/, (NHC6)
/AmmC7/, (NHC7)
/AmmC12/, (NHC12)
/invAb/, (invabasic)
/invdT/, (invdT)
/VPmU/, (vinu)
*, ps
/2-C16U/, 2-C16U 
/2-C18w9U/, 2-C18w9U
/JDi-Palm/, JDi-Palm
/J2-CONC16U/, J2-CONC16U
/J2-C3NC16U/, J2-C3NC16U
/J-C15Ada/, J-C15Ada
/J-2C15AdaU/, J-2C15AdaU
/J-C16NC6/, J-C16NC6
/R2-C6NH-U/, R2-C6NH-U
/J-M1/, J-M1
/J-B1/, J-B1
/J-B2/, J-B2
/J-M2/, J-M2
/2-C16C/, 2-C16C
/2-C16A/, 2-C16A
/2-C16G/, 2-C16G
/R2-C6NH-G/, R2-C6NH-G
/R2-C6NH-C/, R2-C6NH-C
/J2-CONC16A/, J2-CONC16A
/J2-CONC16C/, J2-CONC16C
/J2-CONC16G/, J2-CONC16G
/J2-C15AdaC/, J2-C15AdaC
/J2-M2U/, J2-M2U
/J2-B2U/, J2-B2U
/J2-C3NC16C/, J2-C3NC16C
/J2-C3NC16G/, J2-C3NC16G
/R2-C6NH-A/, R2-C6NH-A
/J2-C15AdaA/, J2-C15AdaA
/J2-C3NC16A/, J2-C3NC16A
/J-C5-SER-1/, J-C5-SER-1
/J-C16-SER-1/, J-C16-SER-1
/J-A2/, J-A2
/J-A1/, J-A1
/J2-C15AdaG/, J2-C15AdaG
/J-C16NAsp/, J-C16NAsp
/J2-C16NC6U/, J2-C16NC6U
/J-C5-REBO-1/, J-C5-REBO-1
/J-C16-REBO-1/, J-C16-REBO-1
/J-C16-IND-1/, J-C16-IND-1
/J-C5-IND-1/, J-C5-IND-1
/J-1C15Ada-2Man/, J-1C15Ada-2Man
/JG-1C15Ada-23DiMan/, JG-1C15Ada-2,3DiMan
/J-TriManPC/, J-TriManPC
/J-triManPO/, J-triManPO
/J-A4/, J-A4
/J-Ara-1/, J-Ara-1
/J-Ara-2/, J-Ara-2
/J-AcCS/, J-AcCS
/J-CbCS/, J-CbCS
/J-MtCD/, J-MtCD`;


const codesWithSmiles = getAllCodesOfSynthesizer(SYNTHESIZERS.GCRS);
const allGcrsCodes = DG.DataFrame.fromCsv(lcmsToGcrs).getCol('GCRS').toList();
export const gcrsCodesWithoutSmiles = differenceOfTwoArrays(allGcrsCodes, codesWithSmiles);
for (const e of gcrsCodesWithoutSmiles)
  map[SYNTHESIZERS.GCRS]['Others'][e] = {name: '', weight: 0, normalized: '', SMILES: ''};


export const weightsObj: {[code: string]: number} = {};
for (const synthesizer of Object.keys(map)) {
  for (const technology of Object.keys(map[synthesizer])) {
    for (const code of Object.keys(map[synthesizer][technology]))
      weightsObj[code] = map[synthesizer][technology][code].weight;
  }
}
for (const [key, value] of Object.entries(MODIFICATIONS))
  weightsObj[key] = value.molecularWeight;
