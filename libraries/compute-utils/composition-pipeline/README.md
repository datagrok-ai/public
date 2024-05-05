# Composition Pipeline

A declarative way to wire multiple script views in a pipeline, as well
as wiring multiple pipelines into a single one.

## Related abstractions overview

- [Function
  call](https://datagrok.ai/help/datagrok/concepts/functions/function-call)
  is constructed from a script or a package function by the platform
  core via [Function
  annotations](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation). It
  provides a runnable object which stores inputs and outputs.

- [RichFunctionView](https://datagrok.ai/help/compute/scripting/Advanced%20scripting/)
  is an abstraction on top of a Function call and is a part of Compute
  package. RichFunctionView is responsible for building UI with inputs
  and outputs widgets and linking them to a function call instance,
  ensuring that once executed, function call inputs are immutable and a
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
cannot be reused directly over several pipelines. Composition pipeline
goal is to provide a way of easily building new pipelines from
existing ones. To achieve this goal, composition pipeline
configuration contains both steps and input/output links
specifications. This allows to compose multiple pipeline
configurations into a single pipeline configuration, which in turn can
be composed with additional pipeline configurations.

## Composition pipeline configuration concepts

- **id** - is an unique name inside a pipeline config.

- **step** - is a description of an individual RichFunctionView
  configuration with some additional customization opportunities.

- **item path** - is a path composed of **ids**. For example, if we
  have a pipeline with id *piplineInner* and a step with id
  *step1*. Inside *piplineInner* **item path** for *step1* is
  `['step1']`. However, if *piplineInner* is used as a nested Pipeline
  to build another pipeline *piplineOuter*, inside this new
  *piplineOuter* the **item path** for *step1* is `['piplineInner',
  'step1']`.

- **link** - is a description of a connection between inputs and
  outputs. It is defined by an array of input **item paths** with the
  corresponding array of output **item paths**. By default input[0]
  value is set to output[0] and so on. Any value change in any of link
  inputs will trigger the link, so an update of all outputs to the
  respective inputs values will happen. Note that links are triggered
  once if multiple inputs are changed synchronously.

- **handler** - is an asynchronous function, if defined it can
  override the default link behavior when a link is triggered. It will
  receive an instance of `RuntimeController` as an argument for
  fetching and setting values. Only a single instance of a link
  handler is active at any point in time, other async handlers running
  will be terminated when using any of RuntimeController methods.

- **validator** - is a link that targets some values as both input and
  output, those inputs are validated. **handler** for such link can
  use `RuntimeController` method `setValidation` and
  `getValidationAction`. Validation data format is the same as for
  [RFV
  validation](https://datagrok.ai/help/compute/scripting/Advanced%20scripting/validating-inputs).

- **popup** - a button in RichFunctionView controls area that launches
  a new arbitrary RichFunctionView. This view is neither a part of
  step sequence nor a part of history, but it can reuse links. Note
  that links from popup will be triggers only when popup is
  confirmed. A path to a popup is prefixed by a step id.

- **action** - a button in RichFunctionView controls area that will
  trigger it's own **handler**, also can be used in **popup** and in
  **validator** actions.

- **hook** - a **handler** that is triggered at a specific pipeline
  livecycle event, like run loading or initialization.

- **compose** is composition pipeline static method, it takes a
  composition configuration and an array of nested pipelines
  configurations. Composition configuration object has additional
  fields to add or remove links, popups and actions for nested
  pipelines, as well as controlling where to put nested pipilines in
  the final pipeline sequence.

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


Now if we want to prefix this pipeline with an another pipeline and
pass its result to input `a` of `AddMock`:


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


A more complicated example will be passing a suggestion for a user to
consider, rather than setting a value directly. In the example below
we suggest to set input `a` of `step2` to the same value as `step1`
`a`. This is is achieved using validators mechanism on `link1`, which
will show a suggestion with a user executable action `action1` which
will set the value.


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
          actions: [{
            id: 'action1',
            friendlyName: 'action1',
            from: ['step1', 'a'],
            to: ['step2', 'a'],
            handler: async ({controller}) => {
              const val = controller.getState(['step1', 'a']);
              controller.setState(['step2', 'a'], val);
            },
            position: 'none',
          }],
        },
      ],
      links: [{
        id: 'link1',
        from: [['step1', 'a'], ['step2', 'a']],
        to: ['step2', 'a'],
        handler: async ({controller}) => {
          const valS = controller.getState(['step1', 'a']);
          const valC = controller.getState(['step2', 'a']);
          if (valC !== valS) {
            const action1 = controller.getValidationAction(['step2', 'action1'], `set to ${valS}`);
            const adv1 = makeAdvice('Try using the provided value', [action1]);
            controller.setValidation(['step2', 'a'], makeValidationResult({notifications: [adv1]}));
          } else
            controller.setValidation(['step2', 'a'], undefined);
        },
      }],
    });
    grok.shell.addView(pipeline.makePipelineView());
    await pipeline.init();
