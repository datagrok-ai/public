/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const DELIMITER = ';'; // what is the need for this?
export const NUCLEOTIDES = ['A', 'G', 'C', 'U', 'T'];

export const TECHNOLOGIES = {
  DNA: 'DNA',
  RNA: 'RNA',
  ASO_GAPMERS: 'ASOGapmers',
  SI_RNA: 'siRNA',
};

export const INPUT_FORMATS = {
  NUCLEOTIDES: 'Nucleotides',
  BIOSPRING: 'BioSpring',
  GCRS: 'GCRS',
  AXOLABS: 'Axolabs',
  MERMADE_12: 'Mermade12',
};

export const SYNTHESIZERS = {
  ...INPUT_FORMATS,
  LCMS: 'LCMS',
};
