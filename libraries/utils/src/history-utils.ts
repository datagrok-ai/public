/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import {DateOptions, FilterOptions} from './function-view';

const getSearchStringByPattern = (datePattern: DateOptions) => {
  switch (datePattern) {
  case 'Today':
    return ` and started > -1d`;
  case 'Yesterday':
    return ` and started > -2d and started < -1d`;
  case 'Any time':
    return ``;
  case 'Last year':
    return `and started > -2y and started < -1y`;
  case 'This year':
    return ` and started > -1y`;
  case 'Last month':
    return ` and started > -2m and started < -1m`;
  case 'This month':
    return ` and started > -1m`;
  case 'Last week':
    return ` and started > -2w and started < -1w`;
  case 'This week':
    return ` and started > -1w`;
  }
};

export namespace historyUtils {
  export async function loadRun(funcCallId: string) {
    const pulledRun = await grok.dapi.functions.calls.allPackageVersions().include('inputs, outputs').find(funcCallId);
    // FIX ME: manually get script since pulledRun contains empty Func
    const script = await grok.dapi.functions.allPackageVersions().find(pulledRun.func.id);
    pulledRun.func = script;
    pulledRun.options['isHistorical'] = true;
    const dfOutputs = wu(pulledRun.outputParams.values() as DG.FuncCallParam[])
      .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const output of dfOutputs)
      pulledRun.outputs[output.name] = await grok.dapi.tables.getTable(pulledRun.outputs[output.name]);

    const dfInputs = wu(pulledRun.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const input of dfInputs)
      pulledRun.inputs[input.name] = await grok.dapi.tables.getTable(pulledRun.inputs[input.name]);

    return pulledRun;
  }

  export async function saveRun(callToSave: DG.FuncCall) {
    const dfOutputs = wu(callToSave.outputParams.values() as DG.FuncCallParam[])
      .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const output of dfOutputs)
      await grok.dapi.tables.uploadDataFrame(callToSave.outputs[output.name]);

    const dfInputs = wu(callToSave.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const input of dfInputs)
      await grok.dapi.tables.uploadDataFrame(callToSave.inputs[input.name]);

    return await grok.dapi.functions.calls.save(callToSave);
  }

  export async function deleteRun(callToDelete: DG.FuncCall) {
    await grok.dapi.functions.calls.delete(callToDelete);
  }

  /**
   * Loads all the function call of this function.
   * Designed to pull hstorical runs in fast manner and the call {@link loadRun} with specified run ID.
   * WARNING: FuncCall inputs/outputs fields are not included
   * @param funcId ID of Func which calls we are looking for. Get it using {@link func.id} field
   * @return Promise on array of FuncCalls corresponding to the passed Func ID
   * @stability Deprecated. Script ID changes with every package release, so searching by ID is useless in practice.
 */
  export async function pullRuns(
    funcId: string,
    filterOptions: FilterOptions = {},
    listOptions: {pageSize?: number, pageNumber?: number, filter?: string, order?: string} = {}
  ): Promise<DG.FuncCall[]> {
    let filteringString = `func.id="${funcId}"`;
    filteringString += filterOptions.author ? ` and session.user.id="${filterOptions.author.id}"`:'';
    filteringString += filterOptions.date ? getSearchStringByPattern(filterOptions.date): '';
    const filter = grok.dapi.functions.calls
      .allPackageVersions()
      .filter(filteringString)
      .include('session.user, options');
    const list = filter.list(listOptions);
    return list;
  }

  /**
   * Loads all the function call of this function.
   * Designed to pull hstorical runs in fast manner and the call {@link loadRun} with specified run ID.
   * WARNING: FuncCall inputs/outputs fields are not included
   * @param funcName Name of Func which calls we are looking for. Get it using {@link func.name} field
   * @return Promise on array of FuncCalls corresponding to the passed Func ID
   * @stability Experimental
 */
  export async function pullRunsByName(
    funcName: string,
    filterOptions: FilterOptions = {},
    listOptions: {pageSize?: number, pageNumber?: number, filter?: string, order?: string} = {}
  ): Promise<DG.FuncCall[]> {
    let filteringString = `func.name="${funcName}"`;
    filteringString += filterOptions.author ? ` and session.user.id="${filterOptions.author.id}"`:'';
    filteringString += filterOptions.date ? getSearchStringByPattern(filterOptions.date): '';
    filteringString += filterOptions.isShared ? ` and options.isShared="true"`: '';
    filteringString += filterOptions.text ?
      ` and ((options.title like "${filterOptions.text}") or (options.annotation like "${filterOptions.text}"))`: '';
    const result =
      grok.dapi.functions.calls
        .allPackageVersions()
        .filter(`(${filteringString})`)
        .include('session.user, options')
        .list(listOptions);
    return result;
  }
}
