// Workflow with validators, meta links, and actions.
// Demonstrates: cross-field validation, UI metadata, and user-triggered actions.
//
// Install: npm i @datagrok-libraries/compute-api

import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';

//name: ValidatedWorkflow
//tags: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function validatedWorkflow(): PipelineConfiguration {
  return {
    id: 'validated',
    nqName: 'MyPackage:ValidatedWorkflow',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'config',
        nqName: 'MyPackage:ConfigureAnalysis',
        friendlyName: 'Configuration',
      },
      {
        id: 'run',
        nqName: 'MyPackage:RunAnalysis',
        friendlyName: 'Run analysis',
      },
    ],
    links: [
      // Data link: pass selected method to the run step
      {
        id: 'method-link',
        from: 'in1:config/method',
        to: 'out1:run/method',
        defaultRestrictions: {out1: 'disabled'},
      },

      // Validator: check that min < max
      {
        id: 'range-validator',
        type: 'validator',
        from: ['in1:run/minValue', 'in2:run/maxValue'],
        to: 'out1:run/minValue',
        handler({controller}) {
          const min = controller.getFirst<number>('in1');
          const max = controller.getFirst<number>('in2');
          if (min != null && max != null && min >= max)
            controller.setValidation('out1', {errors: [{description: 'Min must be less than max'}]});
          else
            controller.setValidation('out1', undefined);
        },
      },

      // Meta link: hide advanced inputs unless expert mode is on
      {
        id: 'expert-mode',
        type: 'meta',
        from: 'in1:config/expertMode',
        to: ['out1:run/advancedParam1', 'out2:run/advancedParam2'],
        handler({controller}) {
          const expert = controller.getFirst<boolean>('in1');
          const meta = expert ? {} : {hidden: true};
          controller.setViewMeta('out1', meta);
          controller.setViewMeta('out2', meta);
        },
      },
    ],

    // Action: reset all inputs to defaults (user-triggered via button)
    actions: [
      {
        id: 'reset',
        from: [],
        to: ['out1:run/minValue', 'out2:run/maxValue'],
        position: 'none',
        friendlyName: 'Reset to defaults',
        handler({controller}) {
          controller.setAll('out1', 0);
          controller.setAll('out2', 100);
        },
      },
    ],
  };
}
