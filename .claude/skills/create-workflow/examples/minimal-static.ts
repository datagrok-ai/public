// Minimal static workflow: two scripts in sequence, no links.
// The user manually fills inputs for each step and runs them.
//
// Install: npm i @datagrok-libraries/compute-api

import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import timezone from 'dayjs/plugin/timezone';

dayjs.extend(utc);
dayjs.extend(timezone);

export const _package = new DG.Package();

//name: MinimalWorkflow
//tags: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function minimalWorkflow(): PipelineConfiguration {
  return {
    id: 'minimal',
    nqName: 'MyPackage:MinimalWorkflow',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'prepare',
        nqName: 'MyPackage:PrepareData',
        friendlyName: 'Prepare data',
      },
      {
        id: 'analyze',
        nqName: 'MyPackage:AnalyzeData',
        friendlyName: 'Analyze',
      },
    ],
  };
}
