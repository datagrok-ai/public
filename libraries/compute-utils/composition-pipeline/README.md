# Composition Pipeline

## Related abstractions overview

- [Function call](https://datagrok.ai/help/datagrok/concepts/functions/function-call)
  is constructed by the platform core via [Function
  annotations](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation). It
  provides a runnable object which stores inputs and outputs.

- [RichFunctionView](https://datagrok.ai/help/compute/scripting-advanced#running-scripts-with-richfunctionview)
  is an abstraction on top of a Function call and is a part of Compute
  package. RichFunctionView is resposible for building UI with inputs
  and outputs widgets and linking them to a function call instance,
  ensuring that once exuted, function call inputs are immutable and a
  fresh funcall is linked to UI inputs. Platform delegates
  scripts/functions inputs filling and running to RichFunctionView
  based on an editor annotation `editor: Compute:RichFunctionViewEditor`.

- __PipelineView__ is a typescript class, located in compute-utils
  library. The main purpose of this class is to manage a sequence of
  RichFunctionView instances (created by Compute package). The
  sequence is defined by constructor arguments and is immutable. Any
  inputs to outputs value passing could be done only by creating a
  PipelineView subclass with custom logic. The resulting UI build by
  PipelineView is a wizard where each step is essentially a
  RichFunctionView.

## Composition pipeline goal

A subclassing based customization is not very flexible and practically
cannot be reused directly over serveral pipelines. Composition
pipeline goal is to provide a way of easily building new piplines from
existing ones. To achive this goal, composition pipeline configuration
contains both steps and input/output links specifications. This allows
to compose multiple pipline configrations in a single pipeline
configuration, which in turn can be composed with additional pipeline
configurations.

## Composition pipeline main concepts

- **id** - is an unique name inside a pipeline config.

- **step** - is a description of an individual RichFunctionView
  configuration with some additional customization opportunities.

- **item path** - is a path composed of **ids**, for example if we
  have a pipeline with id *piplineInner* and a step with id
  *step1*. Inside *piplineInner* **item path** for *step1* is
  `['step1']`. However if *piplineInner* is used as a nested config to
  build another pipeline *piplineOuter*, inside this new
  *piplineOuter* the **item path** for *step1* is
  `['piplineInner', 'step1']`.

- **link** - is a description of a connection between inputs and
  outputs. It is defined by an array of input **item paths** with the
  correspoding array of output **item paths**. By default input[0]
  value is set to output[0] and so on. Any value change in any of link
  inputs will trigger the link, so an update of all outputs to the
  respective inputs values will happen. Note that links are triggered
  once if multiple inputs are changed synchronously.

- **handler** - is an asynchronous function, if defined it can
  override the default link behaviour when a link is triggered. It
  will recieve an instance of `RuntimeController` as an argument for
  fetching and setting values. Only a single instance of a link
  handler is acive at one time, other async handlers running will be
  terminated when using any of RuntimeController methods.

- **popup** - a button in RichFunctionView controlls area that
  launches a new arbitrary RichFunctionView. This view is neither a
  part of step sequence nor a part of history, but it can reuse
  links. Note that links from popup will be triggers only when popup
  is confirmed. A path to a popup is prefixed by a step.

- **action** - a button in RichFunctionView controlls area that will
  triger it's own **handler**, also can be used in **popup**.

- **hook** - a **handler** that is triggered at a specific pipeline
  livecycle event, like run loading or initialization.

- **compose** is composition pipeline static method, it takes a
  composition configuration and an array of nested pipeline
  configurations. Composition configuration object has additional
  fileds to add or remove steps, popups and actions for nestes
  pipelines. The resulting pipeline steps are composition
  configuration steps first then nested pipelines steps.

## Examples

For example, we have two functions in package `MyPackage`:


    //name: AddMock
    //language: javascript
    //input: double a
    //input: double b
    //output: double res
    //editor: Compute:PipelineStepEditor

    res = a + b;



    //name: MulMock
    //language: javascript
    //input: double a
    //input: double b
    //output: double res
    //editor: Compute:PipelineStepEditor

    res = a * b;


If we want to pass result of `AddMock` to argument `a` of `MulMock`, code will
look like:


    const pipeline = new CompositionPipeline({
      id: 'testPipeline',
      nqName: 'MyPackage:MockWrapper1',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:AddMock',
        },
        {
          id: 'step2',
          nqName: 'MyPackage:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['step2', 'a'],
      }]
    });
    grok.shell.addView(pipeline.makePipelineView());
    await pipeline.init();


If we want to pass result of `AddMock` doubled to argument `a` of `MulMock`, code will
look like:


    const pipeline = new CompositionPipeline({
      id: 'testPipeline',
      nqName: 'MyPackage:MockWrapper1',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:AddMock',
        },
        {
          id: 'step2',
          nqName: 'MyPackage:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'a'],
        to: ['step2', 'a'],
        handler: async ({controller}) => {
          const val = controller.getState(['step1', 'a']);
          controller.updateState(['step2', 'a'], val * 2);
        }
      }]
    });
    grok.shell.addView(pipeline.makePipelineView());
    await pipeline.init();


Now we want to prefix this pipeline with a new step and pass its
result to input `a` of `AddMock`:


    //name: MockTripple
    //input: double a
    //output: double res
    export function MockTripple(a: number) {
      return 3 * a;
    }



    const conf1: PipelineConfiguration = {
      id: 'testPipeline1',
      nqName: 'MyPackage:MockWrapper1',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:AddMock',
        },
        {
          id: 'step2',
          nqName: 'MyPackage:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['step2', 'a'],
      }]
    }
    const conf2: PipelineCompositionConfiguration = {
      id: 'testPipeline2',
      nqName: 'MyPackage:MockWrapper10',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:MockTripple',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['testPipeline1', 'step1', 'a'],
      }]

    }
    const conf = CompositionPipeline.compose(conf2, [conf1]);
    const pipeline = new CompositionPipeline(conf);
    grok.shell.addView(pipeline.makePipelineView());
    await pipeline.init();


Now we want to replace `step1` with a new `step1` and pass its result
to `MulMock` input `a`:


    const conf1: PipelineConfiguration = {
      id: 'testPipeline1',
      nqName: 'MyPackage:MockWrapper1',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:AddMock',
        },
        {
          id: 'step2',
          nqName: 'MyPackage:MulMock',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['step2', 'a'],
      }]
    }
    const conf2: PipelineCompositionConfiguration = {
      id: 'testPipeline2',
      nqName: 'MyPackage:MockWrapper11',
      steps: [
        {
          id: 'step1',
          nqName: 'MyPackage:MockTripple',
        },
      ],
      links: [{
        id: 'link1',
        from: ['step1', 'res'],
        to: ['testPipeline1', 'step2', 'a'],
      }],
      itemsToRemove: [['testPipeline1', 'step1'], ['testPipeline1', 'link1']],
    }
    const conf = CompositionPipeline.compose(conf2, [conf1]);
    const pipeline = new CompositionPipeline(conf);
    grok.shell.addView(pipeline.makePipelineView());
    await pipeline.init();
