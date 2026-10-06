// With a JSON schema the result carries `structuredOutput` matching it.

const schema = {
  type: 'object',
  properties: {colours: {type: 'array', items: {type: 'string'}}},
  required: ['colours'],
};

const {structuredOutput} = await grok.ai.prompt('Name three colours.', {schema});
grok.shell.info(structuredOutput.colours.join(', '));
