// A document attachment. Not every engine takes documents; those that cannot throw.

const bytes = await grok.dapi.files.readAsBytes('System:DemoFiles/texts/pdf/research-paper.pdf');
const pdf = new Blob([bytes], {type: 'application/pdf'});
const {text} = await grok.ai.prompt('Summarize this paper in one sentence.',
  {attachments: [{type: 'document', data: pdf, title: 'research-paper.pdf'}]});
grok.shell.info(text);
