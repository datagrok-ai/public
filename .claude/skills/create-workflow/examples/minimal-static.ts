// Minimal static workflow: two scripts in sequence, no links.
// The user manually fills inputs for each step and runs them.

import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';

//name: MinimalWorkflow
//tags: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function minimalWorkflow(): PipelineConfiguration {
  return {
    id: 'minimal',
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
