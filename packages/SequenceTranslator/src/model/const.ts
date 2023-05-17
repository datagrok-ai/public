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

export enum INPUT_FORMATS {
  HELM = 'HELM',
  NUCLEOTIDES = 'Nucleotides',
  BIOSPRING = 'BioSpring',
  GCRS = 'GCRS',
  AXOLABS = 'Axolabs',
  MERMADE_12 = 'Mermade12',
};

export enum SYNTHESIZERS {
  HELM = 'HELM', // helm is not a synthesizer, rename enum
  NUCLEOTIDES = 'Nucleotides',
  BIOSPRING = 'BioSpring',
  GCRS = 'GCRS',
  AXOLABS = 'Axolabs',
  MERMADE_12 = 'Mermade12',
  LCMS = 'LCMS',
};
