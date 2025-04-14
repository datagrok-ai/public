//name: MockProvider3Script
//language: javascript
//output: object config
//editor: Compute2:TreeWizardEditor

config = {
  id: 'pipelinePar',
  nqName: 'Compute2:MockWrapper3', // for history
  provider: 'Compute2:MockProvider3', // for config
  friendlyName: 'Tree wizard model',
  version: '1.0',
  type: 'parallel',
  stepTypes: [{
    type: 'ref',
    provider: 'Compute2:MockProvider2',
    version: '1.0',
  }],
  initialSteps: [
    {
      id: 'pipelinePar',
    },
  ],
};
