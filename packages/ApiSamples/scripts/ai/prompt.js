// One turn in a fresh chat. Every turn ends in an AIResult: text, usage, and structured output
// when a schema was set. Throws when no engine can serve the request.

const {text, usage} = await grok.ai.prompt('Name three colours.');
grok.shell.info(`${text} (${usage.outputTokens ?? '?'} output tokens)`);
