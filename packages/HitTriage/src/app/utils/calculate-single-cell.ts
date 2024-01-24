import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageTemplateFunction, HitTriageTemplateScript} from '../types';
import {HitDesignMolColName, ViDColFormat} from '../consts';
import {joinQueryResults} from '../utils';

export async function calculateSingleCellValues(
  value: string, descriptors: string[], functions: HitTriageTemplateFunction[], scripts: HitTriageTemplateScript[] = [],
  queries: HitTriageTemplateScript[] = [],
): Promise<DG.DataFrame> {
  const canonicalSmiles = grok.chem.convert(value, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
  const col = DG.Column.fromStrings(HitDesignMolColName, [canonicalSmiles]);
  const table = DG.DataFrame.fromColumns([col]);
  await table.meta.detectSemanticTypes();

  if (descriptors.length)
    await grok.chem.descriptors(table, col.name, descriptors);

  for (const func of functions) {
    const props = func.args;
    const fs = DG.Func.find({package: func.package, name: func.name});
    if (!fs.length || !fs[0]) {
      console.warn(`Function ${func.name} from package ${func.package} is not found`);
      continue;
    }
    const f = fs[0];
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

  for (const query of queries) {
    const props = query.args;
    try {
      const queryFunc = await grok.dapi.queries.find(query.id);
      if (queryFunc) {
        const listPropName = queryFunc.inputs[0].name;
        const valueList = [canonicalSmiles];
        const qRes = await queryFunc.apply({...props, [listPropName]: valueList});
        if (!qRes || qRes.rowCount === 0)
          continue;
        await joinQueryResults(table, HitDesignMolColName, qRes);
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
