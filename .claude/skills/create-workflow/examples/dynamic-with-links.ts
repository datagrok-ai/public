// Dynamic workflow with data links and restrictions.
// Users can add/remove analysis steps. The output of
// the prepare step feeds into each analysis step's input.
//
// Install: npm i @datagrok-libraries/compute-api

import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';

//name: DynamicWorkflow
//tags: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function dynamicWorkflow(): PipelineConfiguration {
  return {
    id: 'dynamic-analysis',
    nqName: 'MyPackage:DynamicWorkflow',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'prepare',
        nqName: 'MyPackage:PrepareData',
      },
      {
        id: 'analyses',
        type: 'dynamic',
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
