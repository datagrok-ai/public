import * as DG from 'datagrok-api/dg';
import { MolTrackProp, PROPERTIES, ResultOutput } from './constants';
import { MolTrackProperty } from './types';
import { MolTrackDockerService } from '../services/moltrack-docker-service';

export async function getCorporateCompoundIdByExactStructure(structure: string): Promise<string | null> {
  try {
    const query = {
      'level': 'compounds',
      'output': ['compounds.canonical_smiles', 'compounds.details.corporate_compound_id'],
      'filter': {
        'field': 'compounds.structure',
        'operator': 'IS SIMILAR',
        'value': structure,
        'threshold': 1,
      },
      'output_format': ResultOutput.JSON,
    };
    const { data } = await MolTrackDockerService.search(query, 'compounds');
    return data?.[0]?.['compounds.details.corporate_compound_id'] ?? null;
  } catch (e) {
    console.error('Exact structure search failed:', e);
    return null;
  }
}

export async function getBatchInfoBySynonym(
  synonymName: string, synonymValues: string[], outputFields: string[],
): Promise<any[]> {
  const query = {
    'level': 'compounds',
    'output': outputFields,
    'filter': {
      'field': synonymName,
      'operator': 'IN',
      'value': synonymValues,
    },
    'output_format': ResultOutput.JSON,
  };
  const { data } = await MolTrackDockerService.search(query, 'compounds');
  return data;
}

export function buildPropertyOptions(
  p: MolTrackProp,
  opts?: { reserved?: string[], skipReservedCheck?: boolean },
): any {
  const options: any = {
    name: p.name,
    friendlyName: p.friendly_name,
    type: p.value_type,
  };

  if (p.description)
    options.description = p.description;

  const reserved = opts?.reserved ?? [];
  const skipReservedCheck = opts?.skipReservedCheck ?? false;

  if (p.pattern && (!skipReservedCheck || !reserved.includes(p.name))) {
    const regex = new RegExp(p.pattern);
    options.valueValidators = [
      (val: any) =>
        val == null || val === '' ?
          null :
          regex.test(String(val)) ?
            null :
            `Value does not match pattern: ${p.pattern}`,
    ];
  }

  return options;
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
