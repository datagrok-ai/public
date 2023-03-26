/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const DELIMITER = ';';
export const NUCLEOTIDES = ['A', 'G', 'C', 'U', 'T'];

export const TECHNOLOGIES = {
  DNA: 'DNA',
  RNA: 'RNA',
  ASO_GAPMERS: 'For ASO Gapmers',
  SI_RNA: 'For 2\'-OMe and 2\'-F modified siRNA',
};

export const INPUT_FORMATS = {
  RAW_NUCLEOTIDES: 'Raw Nucleotides',
  BIOSPRING: 'BioSpring Codes',
  GCRS: 'Janssen GCRS Codes',
  AXOLABS: 'Axolabs Codes',
  MERMADE_12: 'Mermade 12',
};

export const SYNTHESIZERS = {
  ...INPUT_FORMATS,
  LCMS: 'LCMS',
};
