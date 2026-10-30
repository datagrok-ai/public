/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TemplateFunction, TemplateScript, IComputeDialogResult} from './types';
import {funcTypeNames} from './consts';
import {joinQueryResults} from './utils';

export type ComputeExecutionOptions = {
  convertToSmiles?: (value: string) => string,
  logger?: {error: (...args: any[]) => void},
  molColName?: string,
}

export async function calculateCellValues(
  values: string[], descriptors: string[], functions: TemplateFunction[], scripts: TemplateScript[] = [],
  queries: TemplateScript[] = [],
  options?: ComputeExecutionOptions,
): Promise<DG.DataFrame> {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  const convertToSmiles = options?.convertToSmiles ?? ((v: string) => v);
  const logger = options?.logger ?? console;
  const molColName = options?.molColName ?? 'Molecule';
  const canonicalSmiles = values.map((v) => convertToSmiles(v));
  const col = DG.Column.fromStrings(molColName, canonicalSmiles);
  const table = DG.DataFrame.fromColumns([col]);
  table.name = 'HD cell values';
  await table.meta.detectSemanticTypes();
  const promises: Promise<void>[] = [];
  promises.push(new Promise(async (resolve) => {
    try {
      if (descriptors.length)
        await grok.chem.descriptors(table, col.name, descriptors);
    } catch (e) {
      logger.error('Descriptors calculation error');
      logger.error(e);
    }
    resolve();
  }));

  for (const func of functions) {
    promises.push(new Promise(async (resolve) => {
      try {
        const props = func.args;
        const fs = DG.Func.find({package: func.package, name: func.name});
        if (!fs.length || !fs[0]) {
          console.warn(`Function ${func.name} from package ${func.package} is not found`);
          resolve();
          return;
        }
        const f = fs[0];
        const tablePropName = f.inputs[0].name;
        const colPropName = f.inputs[1].name;
        await f.apply({...props, [tablePropName]: table, [colPropName]: col.name});
      } catch (e) {
        logger.error(e);
      }
      resolve();
    }));
  }

  for (const script of scripts) {
    const props = script.args;
    promises.push(new Promise(async (resolve) => {
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
        logger.error(e);
      }
      resolve();
    }));
  }

  for (const query of queries) {
    const props = query.args;
    promises.push(new Promise(async (resolve) => {
      try {
        const queryFunc = DG.Func.find({name: query.name})
          .filter((s) => s.type === funcTypeNames.query && s.id === query.id)[0];
        if (queryFunc) {
          const listPropName = queryFunc.inputs[0].name;
          const valueList = [canonicalSmiles];
          const qRes = await queryFunc.apply({...props, [listPropName]: valueList});
          if (!qRes || qRes.rowCount === 0) {
            resolve();
            return;
          }
          await joinQueryResults(table, molColName, qRes);
        }
      } catch (e) {
        logger.error(e);
      }
      resolve();
    }));
  }
  await Promise.all(promises);
  pg.close();
  const molCol = table.col(molColName)!;
  for (let i = 0; i < molCol.length; i++) {
    // restore previous values
    molCol.set(i, values[i], false);
  }
  return table;
}

/** Runs on a copy of the table's listed columns (none hidden behind `~`), the molecule column holding the canonical SMILES under its own name, as `calculateCellValues` does; what the functions add or change fills the columns of the same name, and the table says its values changed. */
export async function calculateColumns(resultMap: IComputeDialogResult, dataFrame: DG.DataFrame, molColName: string, options?: ComputeExecutionOptions) {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  const convertToSmiles = options?.convertToSmiles ?? ((v: string) => v);
  const logger = options?.logger ?? console;
  const molCol = dataFrame.col(molColName);
  if (!molCol)
    throw new Error('There is no molecule column in dataframe');

  const molSmilesCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, molColName,
    molCol.toList().map((value, i) => molCol.isNone(i) || !value ? null : convertToSmiles(value)));
  molSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  const table = DG.DataFrame.fromColumns([molSmilesCol, ...dataFrame.columns.toList()
    .filter((c) => c.name && !c.name.startsWith('~') && c.name !== molColName).map((c) => c.clone())]);
  const given = new Map(table.columns.toList().map((c) => [c.name, {dart: c.dart, version: c.version}]));
  // run again, a function names its column anew beside the one it made before: "Acceptor (Molecule) (2)" fills "Acceptor (Molecule)"
  const renamed = new Map([...given.keys()].filter((n) => n !== molColName).map((n) => [table.columns.getUnusedName(n), n]));
  const promises: Promise<void>[] = [];
  promises.push(new Promise(async (resolve) => {
    try {
      if (resultMap.descriptors && resultMap.descriptors.length > 0)
        await grok.chem.descriptors(table, molColName, resultMap.descriptors);
    } catch (e) {
      logger.error('Descriptors calculation error');
      logger.error(e);
    }
    resolve();
  }));

  for (const funcName of Object.keys(resultMap.externals)) {
    promises.push(new Promise(async (resolve) => {
      try {
        const props = resultMap.externals[funcName];
        const f = DG.Func.find({package: funcName.split(':')[0], name: funcName.split(':')[1]})[0];
        const tablePropName = f.inputs[0].name;
        const colPropName = f.inputs[1].name;
        if (props)
          await f.apply({...props, [tablePropName]: table, [colPropName]: molColName});
      } catch (e) {
        logger.error(e);
      }
      resolve();
    }));
  };
  // handling scripts
  for (const scriptName of Object.keys(resultMap.scripts ?? {})) {
    const props = resultMap.scripts![scriptName];
    promises.push(new Promise(async (resolve) => {
      if (props) {
        const scriptParts = scriptName.split(':');
        const scriptId = scriptParts[2];
        if (!scriptId) {
          resolve();
          return;
        }
        try {
          const sn = scriptParts[1]; // script name
          const s = DG.Func.find({name: sn})
            .filter((s) => s.type === funcTypeNames.script && s.id === scriptId)[0] as DG.Script | undefined;
          if (!s) {
            resolve();
            return;
          }
          const tablePropName = s!.inputs[0].name;
          const colPropName = s!.inputs[1].name;
          const r: DG.DataFrame = await s.apply({...props, [tablePropName]: table, [colPropName]: molColName});
          if (r && r.rowCount === table.rowCount && s.language === 'python') {
            for (const c of r.columns) {
              c.name = table.columns.getUnusedName(c.name);
              table.columns.add(c);
            }
          }
        } catch (e) {
          logger.error(e);
        }
      }
      resolve();
    }));
  };

  // handling queries
  for (const queryName of Object.keys(resultMap.queries ?? {})) {
    const props = resultMap.queries![queryName];
    promises.push(new Promise(async (resolve) => {
      if (props) {
        const queryParts = queryName.split(':');
        const queryId = queryParts[2];
        if (!queryId) {
          resolve();
          return;
        }
        try {
          const qn = queryParts[1]; // query name
          const s = DG.Func.find({name: qn})
            .filter((s) => s.type === funcTypeNames.query && s.id === queryId)[0];
          if (!s) {
            resolve();
            return;
          }
          const listPropName = s.inputs[0].name;
          const resDf: DG.DataFrame = await s.apply({...props, [listPropName]: molSmilesCol.toList()});
          if (resDf)
            await joinQueryResults(table, molColName, resDf);
        } catch (e) {
          console.error(e);
        }
      }
      resolve();
    }));
  };
  await Promise.all(promises);
  pg.close();
  for (const col of table.columns.toList()) {
    const was = given.get(col.name);
    if (col.name === molColName || (was !== undefined && was.dart === col.dart && was.version === col.version))
      continue;
    const name = was ? col.name : renamed.get(col.name) ?? col.name;
    const target = dataFrame.col(name);
    // a number column takes counts in place: replacing it would take what is bound to it (a filter) away
    if (target?.type === col.type || (target?.type === DG.COLUMN_TYPE.FLOAT && col.type === DG.COLUMN_TYPE.INT)) {
      for (let i = 0; i < col.length; i++)
        target.set(i, col.isNone(i) ? null : col.get(i), false);
      if (col.semType)
        target.semType = col.semType;
      continue;
    }
    const added = col.clone();
    added.name = name;
    if (!target) {
      dataFrame.columns.add(added);
      continue;
    }
    // a column of another type is a new column, which the platform's viewers and filters let go of: they get it back
    const bound = settingsNaming(dataFrame, name);
    dataFrame.columns.replace(target, added);
    for (const [viewer, settings] of bound)
      viewer.setOptions(settings);
  }
  dataFrame.fireValuesChanged();
}

// a Filters viewer takes its `filters` as the whole list, so the list is kept whole
function settingsNaming(dataFrame: DG.DataFrame, name: string): [DG.Viewer, {[key: string]: unknown}][] {
  const found: [DG.Viewer, {[key: string]: unknown}][] = [];
  for (const view of grok.shell.tableViews) {
    if (view.dataFrame?.dart !== dataFrame.dart)
      continue;
    for (const viewer of view.viewers) {
      const settings = Object.fromEntries(Object.entries(viewer.getOptions(true).look ?? {}).filter(([key, value]) =>
        key.endsWith('ColumnName') ? value === name :
          key.endsWith('ColumnNames') ? Array.isArray(value) && value.includes(name) :
            key === 'filters' && Array.isArray(value) && value.some((f) => f?.column === name)));
      if (Object.keys(settings).length > 0)
        found.push([viewer, settings]);
    }
  }
  return found;
}
