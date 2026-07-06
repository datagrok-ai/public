// Flows (visual workflows from the Flow package) are Script entities with language 'flow':
// the body is an annotation header followed by the flow's JSON document.
// Saving one creates a regular platform entity — shareable, searchable, runnable.
(async () => {
  if (DG.Func.find({package: 'Flow', name: 'flowScriptHandler'}).length === 0) {
    grok.shell.info('Install the Flow package to run this sample');
    return;
  }

  const name = 'DemoFlow_' + Math.random().toString(36).substring(7);
  // The smallest runnable flow: an int input passed through to an int output.
  const body = `//name: ${name}
//language: flow
//tags: flow
//input: int a = 3
//output: int result
{"version":"2.0","name":"${name}","description":"","author":"","created":"","modified":"",
"nodes":[
  {"id":"in","typeName":"Inputs/Int Input","label":"Int Input","pos":{"x":0,"y":0},
   "properties":{"paramName":"a","defaultValue":3},"inputValues":{}},
  {"id":"out","typeName":"Outputs/Value Output","label":"Value Output","pos":{"x":300,"y":0},
   "properties":{"paramName":"result","outputType":"int"},"inputValues":{}}],
"connections":[{"id":"c1","source":"in","sourceOutput":"value","target":"out","targetInput":"value"}],
"metadata":{"settings":{"scriptName":"${name}","scriptDescription":"","tags":["flow"]}}}`;

  let flow = DG.Script.create(body);
  flow = await grok.dapi.scripts.save(flow);
  console.log(`Saved flow entity: ${flow.nqName} (language: ${flow.language})`);

  // Runs like any function: the Flow package compiles the graph and executes it.
  const result = await grok.functions.call(flow.nqName, {a: 5});
  console.log(`Flow result: ${result}`);

  // Flows have their own gallery (Browse | Platform | Functions | Flows) and
  // open in the visual editor on double-click.
  await grok.dapi.scripts.delete(flow);
  grok.shell.info(`Flow demo completed: ${name} returned ${result}`);
})();
