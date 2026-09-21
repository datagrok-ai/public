// A document attachment. Not every engine takes documents; those that cannot throw.

const text = await grok.dapi.files.readAsText('System:DemoFiles/texts/python.txt');
const doc = new Blob([text], {type: 'text/plain'});
const {text: summary} = await grok.ai.prompt('Summarize this document in one sentence.',
  {engine: 'gemma', attachments: [{type: 'document', data: doc, title: 'python.txt'}]});
grok.shell.info(summary);
