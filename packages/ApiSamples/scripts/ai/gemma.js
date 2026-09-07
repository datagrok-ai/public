// Gemma through Chrome's built-in model (Chrome 138+, registered by Grokky as the 'gemma' engine).
// Everything runs on this machine: the model is downloaded once, then works offline at no per-token cost.
// Browsers without it report the engine as unavailable, and calls that do not name an engine skip it.

const gemma = await grok.ai.getEngine('gemma');
if (!gemma)
  throw new Error('The gemma engine is not registered: install Grokky');
if (!(await gemma.available()))
  await gemma.download();

const pause = () => DG.delay(5000);
const opts = {engine: 'gemma'};
const subs = [
  grok.ai.onPrompt.subscribe((e) => grok.shell.info(`prompt: ${e.prompt}`)),
  grok.ai.onResult.subscribe((e) => grok.shell.info(`result: ${e.result.text.length} chars`)),
  grok.ai.onCancel.subscribe(() => grok.shell.info('cancelled')),
  grok.ai.onError.subscribe((e) => grok.shell.warning(e.error.message)),
];
try {
  // One turn
  grok.shell.info((await grok.ai.prompt('Name three colours.', opts)).text);
  await pause();

  // Streaming, cancelled from the AbortSignal after the first few deltas
  const ac = new AbortController();
  let streamed = '';
  try {
    for await (const chunk of grok.ai.stream('Write a long story about a dog.', {...opts, signal: ac.signal})) {
      streamed += chunk;
      if (streamed.length > 80)
        ac.abort();
    }
  } catch (e) {
    if (!ac.signal.aborted)
      throw e;
  }
  grok.shell.info(`streamed before cancel: ${streamed}`);
  await pause();

  // Structured output: responseConstraint under the hood, so the JSON matches the schema
  const schema = {
    type: 'object',
    properties: {colours: {type: 'array', items: {type: 'string'}}},
    required: ['colours'],
  };
  const {structuredOutput} = await grok.ai.prompt('Name three colours.', {...opts, schema});
  grok.shell.info(structuredOutput.colours.join(', '));
  await pause();

  // An image goes to the model as is
  const bytes = await grok.dapi.files.readAsBytes('System:DemoFiles/images/cats/dog1.jpg');
  const image = new Blob([bytes], {type: 'image/jpeg'});
  grok.shell.info((await grok.ai.prompt('What animal is this? One word.',
    {...opts, attachments: [{type: 'image', data: image}]})).text);
  await pause();

  // A text document is inlined into the prompt; binary documents (PDF) throw with this engine
  const text = await grok.dapi.files.readAsText('System:DemoFiles/texts/python.txt');
  const doc = new Blob([text], {type: 'text/plain'});
  grok.shell.info((await grok.ai.prompt('Summarize this document in one sentence.',
    {...opts, attachments: [{type: 'document', data: doc, title: 'python.txt'}]})).text);
  await pause();

  // Multi-turn: the session keeps the transcript
  const chat = await grok.ai.chat({...opts, systemPrompt: 'Answer in one sentence.'});
  try {
    grok.shell.info((await chat.prompt('Which dose-response curve shape indicates a problem?')).text);
    grok.shell.info((await chat.prompt('Why that one in particular?')).text);
  } finally {
    chat.close();
  }
} finally {
  subs.forEach((s) => s.unsubscribe());
}
