import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISubstructProvider} from '@datagrok-libraries/chem-meta/src/types';

import {ISeqMonomer} from '../helm/types';

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
