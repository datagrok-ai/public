//name: MockPipeline3
//language: javascript
//input: object params
//output: object config
//editor: Compute2:TreeWizardEditor
//tags: test, compute2

config = {
  id: 'pipelinePar',
  nqName: 'Compute2:MockPipeline3',
  friendlyName: 'Pipeline 3',
  version: '1.0',
  type: 'parallel',
  stepTypes: [{
    type: 'ref',
    provider: 'Compute2:MockPipeline2',
    version: '1.0',
  }],
  initialSteps: [
    {
      id: 'pipelinePar',
    },
  ],
};
