/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {
  calculateCellValues as _calculateCellValues,
  calculateColumns as _calculateColumns,
} from '@datagrok-libraries/statistics/src/compute-functions/execution';
import type {ComputeExecutionOptions} from '@datagrok-libraries/statistics/src/compute-functions/execution';
import {TemplateFunction, TemplateScript, IComputeDialogResult} from '../types';
import {HitDesignMolColName} from '../consts';
import {_package} from '../../package';

function getOptions(): ComputeExecutionOptions {
  return {
    convertToSmiles: (v: string) => _package.convertToSmiles(v),
    logger: _package.logger,
    molColName: HitDesignMolColName,
  };
}

export async function calculateCellValues(
  values: string[], descriptors: string[], functions: TemplateFunction[], scripts: TemplateScript[] = [],
  queries: TemplateScript[] = [],
): Promise<DG.DataFrame> {
  return _calculateCellValues(values, descriptors, functions, scripts, queries, getOptions());
}

export async function calculateColumns(resultMap: IComputeDialogResult, dataFrame: DG.DataFrame, molColName: string) {
  return _calculateColumns(resultMap, dataFrame, molColName, getOptions());
}
