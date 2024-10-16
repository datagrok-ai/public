/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const NUCLEOTIDES = ['A', 'G', 'C', 'U'];

export const TECHNOLOGIES = {
  DNA: 'DNA',
  RNA: 'RNA',
  ASO_GAPMERS: 'ASOGapmers',
  SI_RNA: 'siRNA',
};

export enum DEFAULT_FORMATS {
  HELM = 'HELM',
  AXOLABS = 'Axolabs',
}
