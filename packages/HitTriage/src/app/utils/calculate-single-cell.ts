import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageTemplateFunction, HitTriageTemplateScript, IComputeDialogResult} from '../types';
import {funcTypeNames, HitDesignMolColName, ViDColFormat} from '../consts';
import {joinQueryResults} from '../utils';

export async function calculateSingleCellValues(
  value: string, descriptors: string[], functions: HitTriageTemplateFunction[], scripts: HitTriageTemplateScript[] = [],
  queries: HitTriageTemplateScript[] = [],
): Promise<DG.DataFrame> {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  // TODO: this converts value to canonical one. We need to do it in a better way
  const canonicalSmiles = grok.chem.convert(value, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
  const col = DG.Column.fromStrings(HitDesignMolColName, [canonicalSmiles]);
  const table = DG.DataFrame.fromColumns([col]);
  table.name = 'HD Single cell values';
  await table.meta.detectSemanticTypes();

  try {
    if (descriptors.length)
      await grok.chem.descriptors(table, col.name, descriptors);
  } catch (e) {
    console.error('Descriptors calculation error', e);
  }
  for (const func of functions) {
    try {
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
    } catch (e) {
      console.error(e);
      continue;
    }
  }

  for (const script of scripts) {
    const props = script.args;
    try {
      const scriptFunc = DG.Func.find({name: script.name})
        .filter((s) => s.type === funcTypeNames.script && s.id === script.id)[0] as DG.Script | undefined;
      if (scriptFunc) {
        const tablePropName = scriptFunc.inputs[0].name;
        const colPropName = scriptFunc.inputs[1].name;

        const r: DG.DataFrame = await scriptFunc.apply({...props, [tablePropName]: table, [colPropName]: col.name});
        if (r && r.rowCount === table.rowCount && scriptFunc.language === 'python') {
          for (const c of r.columns) {
            c.name = table.columns.getUnusedName(c.name);
            table.columns.add(c);
          }
        }
      }
    } catch (e) {
      console.error(e);
    }
  }

  for (const query of queries) {
    const props = query.args;
    try {
      const queryFunc = DG.Func.find({name: query.name})
        .filter((s) => s.type === funcTypeNames.query && s.id === query.id)[0];
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
  pg.close();
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

export async function calculateColumns(resultMap: IComputeDialogResult, dataFrame: DG.DataFrame, molColName: string) {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  // first step: convert all values to canonical smiles.
  const molCol = dataFrame.col(molColName);
  if (!molCol)
    throw new Error('There is no molecule column in dataframe');
  for (let i = 0; i < molCol.length; i++) {
    if (molCol.isNone(i))
      continue;
    const value: string | null = molCol.get(i);
    if (!value)
      continue;
    const newVal = grok.chem.convert(value, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
    molCol.set(i, newVal, false);
  }
  try {
    if (resultMap.descriptors && resultMap.descriptors.length > 0)
      await grok.chem.descriptors(dataFrame!, molColName!, resultMap.descriptors);
  } catch (e) {
    console.error('Descriptors calculation error', e);
  }

  for (const funcName of Object.keys(resultMap.externals)) {
    try {
      const props = resultMap.externals[funcName];
      const f = DG.Func.find({package: funcName.split(':')[0], name: funcName.split(':')[1]})[0];
      const tablePropName = f.inputs[0].name;
      const colPropName = f.inputs[1].name;
      if (props)
        await f.apply({...props, [tablePropName]: dataFrame!, [colPropName]: molColName});
    } catch (e) {
      console.error(e);
    }
  };
  // handling scripts
  for (const scriptName of Object.keys(resultMap.scripts ?? {})) {
    const props = resultMap.scripts![scriptName];
    if (props) {
      const scriptParts = scriptName.split(':');
      const scriptId = scriptParts[2];
      if (!scriptId)
        continue;
      try {
        const sn = scriptParts[1]; // script name
        const s = DG.Func.find({name: sn})
          .filter((s) => s.type === funcTypeNames.script && s.id === scriptId)[0] as DG.Script | undefined;
        if (!s)
          continue;
        const tablePropName = s.inputs[0].name;
        const colPropName = s.inputs[1].name;
        const r: DG.DataFrame = await s.apply({...props, [tablePropName]: dataFrame, [colPropName]: molColName});
        if (r && r.rowCount === dataFrame.rowCount && s.language === 'python') {
          for (const c of r.columns) {
            c.name = dataFrame.columns.getUnusedName(c.name);
            dataFrame.columns.add(c);
          }
        }
      } catch (e) {
        console.error(e);
      }
    }
  };

  // handling queries
  for (const queryName of Object.keys(resultMap.queries ?? {})) {
    const props = resultMap.queries![queryName];
    if (props) {
      const queryParts = queryName.split(':');
      const queryId = queryParts[2];
      if (!queryId)
        continue;
      try {
        const qn = queryParts[1]; // query name
        const s = DG.Func.find({name: qn})
          .filter((s) => s.type === funcTypeNames.query && s.id === queryId)[0];
        if (!s)
          continue;
        const listPropName = s.inputs[0].name;
        const molList = dataFrame!.col(molColName!)!.toList();
        const resDf: DG.DataFrame = await s.apply({...props, [listPropName]: molList});
        if (resDf)
          await joinQueryResults(dataFrame!, molColName!, resDf);
      } catch (e) {
        console.error(e);
      }
    }
  };
  pg.close();
}
