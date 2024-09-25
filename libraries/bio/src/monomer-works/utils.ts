import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISubstructProvider} from '@datagrok-libraries/chem-meta/src/types';
import {PolymerType} from '@datagrok-libraries/js-draw-lite/src/types/org';

import {OrgType} from '../helm/types';
import {ISeqMonomer} from '../helm/types';
import {HELM_POLYMER_TYPE} from '../utils/const';
import {ALPHABET} from '../utils/macromolecule';

declare const org: OrgType;

export type MonomerHoverLink = {
  targetCol: DG.Column,
  handler(seqGridCell: DG.GridCell, monomer: ISeqMonomer | null, targetGridCol: DG.GridColumn): boolean;
} & ISubstructProvider;

export function getUnusedColName(df: DG.DataFrame | undefined, colName: string): string {
  if (!df) return colName;
  return df.columns.getUnusedName(colName);
}

export function getMolColName(df: DG.DataFrame | undefined, seqColName: string): string {
  return getUnusedColName(df, `molfile(${seqColName})`);
}

export function alphabetToPolymerType(alphabet: ALPHABET): PolymerType {
  // determine the polymer type according to HELM specifications
  let resPolymerType: HELM_POLYMER_TYPE;
  // todo: an exception from dart comes before this check if the alphabet is UN
  if (alphabet === ALPHABET.PT || alphabet === ALPHABET.UN)
    resPolymerType = HELM_POLYMER_TYPE.PEPTIDE;
  else if (alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA)
    resPolymerType = HELM_POLYMER_TYPE.RNA;
  else {
    throw new Error(`Unexpected alphabet '${alphabet}'.`);
  }
  return resPolymerType;
}

/** @deprecated Use DG.Color.hexToPercentRgb */
export function hexToPercentRgb(hex: string): number[] | null {
  const result = hex.length === 7 ? /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex) :
    /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result ? [
    parseInt(result[1], 16) / 256,
    parseInt(result[2], 16) / 256,
    parseInt(result[3], 16) / 256,
    result.length > 4 ? parseInt(result[4], 16) / 256 : 0.3
  ] : null;
}
