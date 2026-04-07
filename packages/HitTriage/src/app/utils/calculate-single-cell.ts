/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageTemplateFunction, HitTriageTemplateScript, IComputeDialogResult} from '../types';
import {funcTypeNames, HitDesignMolColName} from '../consts';
import {joinQueryResults} from '../utils';
import {_package} from '../../package';

export async function calculateCellValues(
  values: string[], descriptors: string[], functions: HitTriageTemplateFunction[], scripts: HitTriageTemplateScript[] = [],
  queries: HitTriageTemplateScript[] = [],
): Promise<DG.DataFrame> {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  // TODO: this converts value to canonical one. We need to do it in a better way
  const canonicalSmiles = values.map((v) => _package.convertToSmiles(v));
  const col = DG.Column.fromStrings(HitDesignMolColName, canonicalSmiles);
  const table = DG.DataFrame.fromColumns([col]);
  table.name = 'HD cell values';
  await table.meta.detectSemanticTypes();
  const promises: Promise<void>[] = [];
  promises.push(new Promise(async (resolve) => {
    try {
      if (descriptors.length)
        await grok.chem.descriptors(table, col.name, descriptors);
    } catch (e) {
      _package.logger.error('Descriptors calculation error');
      _package.logger.error(e);
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
        _package.logger.error(e);
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
        _package.logger.error(e);
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
          await joinQueryResults(table, HitDesignMolColName, qRes);
        }
      } catch (e) {
        _package.logger.error(e);
      }
      resolve();
    }));
  }
  await Promise.all(promises);
  pg.close();
  const molCol = table.col(HitDesignMolColName)!;
  for (let i = 0; i < molCol.length; i++) {
    // restore previous values
    molCol.set(i, values[i], false);
  }
  return table;
};

export async function calculateColumns(resultMap: IComputeDialogResult, dataFrame: DG.DataFrame, molColName: string) {
  const pg = DG.TaskBarProgressIndicator.create('Calculating ...');
  // first step: convert all values to canonical smiles.
  const molCol = dataFrame.col(molColName);
  if (!molCol)
    throw new Error('There is no molecule column in dataframe');

  const molSmilesColName = `~${molColName}.canonicalSmiles`; // chem package might have already created that, but we update it just in any case
  const molSmilesCol = dataFrame.columns.getOrCreate(molSmilesColName, DG.TYPE.STRING);
  // we don't want to overwrite existing values, as they can contain important conformation information
  molSmilesCol.semType = DG.SEMTYPE.MOLECULE;

  const molColList = molCol.toList();
  for (let i = 0; i < molColList.length; i++) {
    if (molCol.isNone(i))
      continue;
    const value: string | null = molColList[i];
    if (!value)
      continue;
    const newVal = _package.convertToSmiles(value);
    molSmilesCol.set(i, newVal, false);
  }
  const promises: Promise<void>[] = [];
  promises.push(new Promise(async (resolve) => {
    try {
      if (resultMap.descriptors && resultMap.descriptors.length > 0)
        await grok.chem.descriptors(dataFrame!, molSmilesColName!, resultMap.descriptors);
    } catch (e) {
      _package.logger.error('Descriptors calculation error');
      _package.logger.error(e);
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
          await f.apply({...props, [tablePropName]: dataFrame!, [colPropName]: molSmilesColName});
      } catch (e) {
        _package.logger.error(e);
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
          const r: DG.DataFrame = await s.apply({...props, [tablePropName]: dataFrame, [colPropName]: molSmilesColName});
          if (r && r.rowCount === dataFrame.rowCount && s.language === 'python') {
            for (const c of r.columns) {
              c.name = dataFrame.columns.getUnusedName(c.name);
              dataFrame.columns.add(c);
            }
          }
        } catch (e) {
          _package.logger.error(e);
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
            await joinQueryResults(dataFrame!, molSmilesColName!, resDf);
        } catch (e) {
          console.error(e);
        }
      }
      resolve();
    }));
  };
  await Promise.all(promises);
  pg.close();
}
