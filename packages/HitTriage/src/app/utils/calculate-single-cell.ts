import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageTemplateFunction, HitTriageTemplateScript} from '../types';
import {HitDesignMolColName, ViDColFormat} from '../consts';

export async function calculateSingleCellValues(
  value: string, descriptors: string[], functions: HitTriageTemplateFunction[], scripts: HitTriageTemplateScript[] = [],
): Promise<DG.DataFrame> {
  const col = DG.Column.fromStrings(HitDesignMolColName, [value]);
  const table = DG.DataFrame.fromColumns([col]);
  await table.meta.detectSemanticTypes();

  if (descriptors.length)
    await grok.chem.descriptors(table, col.name, descriptors);

  for (const func of functions) {
    const props = func.args;
    const f = await grok.functions.find(`${func.package}:${func.name}`);
    if (!f)
      continue;
    const tablePropName = f.inputs[0].name;
    const colPropName = f.inputs[1].name;
    await f.apply({...props, [tablePropName]: table, [colPropName]: col.name});
  }

  for (const script of scripts) {
    const props = script.args;
    try {
      const scriptFunc = await grok.dapi.scripts.find(script.id);
      if (scriptFunc) {
        const tablePropName = scriptFunc.inputs[0].name;
        const colPropName = scriptFunc.inputs[1].name;
        await scriptFunc.apply({...props, [tablePropName]: table, [colPropName]: col.name});
      }
    } catch (e) {
      console.error(e);
    }
  }
  return table;
};

export function getNewVid(vidCol: DG.Column<any>) {
  let maxId = 0;
  if (vidCol.length > 1) {
    for (const vid of vidCol.toList()) {
      if (!vid || !vid.startsWith(ViDColFormat[0]) || vid.length !== ViDColFormat.length)
        continue;
      const num = parseInt(vid.substring(1));
      if (isNaN(num))
        continue;
      maxId = Math.max(num, maxId);
    };
  }

  return `V${''.concat(...new Array(ViDColFormat.length - 1 - (maxId +1).toString().length).fill('0'))}${maxId + 1}`;
}
