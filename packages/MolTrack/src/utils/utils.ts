import * as DG from 'datagrok-api/dg';
import { PROPERTIES } from './constants';
import { MolTrackProperty } from './types';

export async function checkCompoundExists(smiles: string): Promise<boolean> {
  return false;
}

export function flattened(item: any, props: DG.Property[]) {
  const row: any = {};
  for (const [key, value] of Object.entries(item)) {
    if (typeof value === 'object' && value !== null) {
      //handle properties object
      if (Array.isArray(value) && key === PROPERTIES) {
        props.forEach((prop: DG.Property) => {
          const valIdx = value.findIndex((it) => it.name.toLowerCase() === prop.name.toLocaleLowerCase());
          const val = valIdx === -1 ? null :
            (value[valIdx] as MolTrackProperty).value_num ??
          (value[valIdx] as MolTrackProperty).value_datetime ??
          (value[valIdx] as MolTrackProperty).value_uuid ??
          (value[valIdx] as MolTrackProperty).value_string;
          row[prop.friendlyName ?? prop.name] = val;
        });
      } else
        row[key] = JSON.stringify(value);
    } else
      row[key] = value;
  }
  return row;
}

export function safeParse(value: any, fallback = []) {
  if (typeof value === 'string') {
    try {
      return JSON.parse(value);
    } catch {
      return fallback;
    }
  }
  return value || fallback;
}
