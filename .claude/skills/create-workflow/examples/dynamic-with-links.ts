// Dynamic workflow with data links and restrictions.
// Users can add/remove analysis steps. The output of
// the prepare step feeds into each analysis step's input.

import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';

//name: DynamicWorkflow
//tags: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function dynamicWorkflow(): PipelineConfiguration {
  return {
    id: 'dynamic-analysis',
    type: 'static',
    steps: [
      {
        id: 'prepare',
        nqName: 'MyPackage:PrepareData',
      },
      {
        id: 'analyses',
        type: 'dynamic',    // 'parallel' and 'sequential' are also accepted
        friendlyName: 'Analyses',
        stepTypes: [
          {
            id: 'regression',
            nqName: 'MyPackage:RunRegression',
            friendlyName: 'Regression',
          },
          {
            id: 'clustering',
            nqName: 'MyPackage:RunClustering',
            friendlyName: 'Clustering',
          },
        ],
        initialSteps: [
          {id: 'regression'},
        ],
      },
    ],
    links: [
      {
        // Pass prepared data to each analysis step's input.
        // all() matches every instance of the step type.
        id: 'data-link',
        from: 'in1:prepare/result',
        to: 'out1:analyses/all(regression|clustering)/data',
        defaultRestrictions: {
          out1: 'restricted',   // user can edit but sees a warning
        },
      },
      {
        // Custom handler: transform the value before passing it
        id: 'params-link',
        from: 'in1:prepare/config',
        to: 'out1:analyses/all(regression|clustering)/params',
        handler({controller}) {
          const config = controller.getFirst('in1');
          // Extract relevant params for each analysis type
          controller.setAll('out1', config?.analysisParams, 'disabled');
        },
      },
    ],
  };
}
