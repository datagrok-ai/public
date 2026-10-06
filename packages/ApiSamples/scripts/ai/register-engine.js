//help-url: https://datagrok.ai/help/develop/function-roles
// A session-scoped engine (packages publish one via `meta.role: aiEngine`).
// An engine extends DG.AIEngine; its chat extends DG.AIChat and implements respond():
// an async generator that yields text deltas and returns the AIResult.

class EchoChat extends DG.AIChat {
  constructor(options) {
    super();
    this.transcript = options.systemPrompt ? [options.systemPrompt] : [];
  }

  async* respond(prompt, options) {
    this.transcript.push(prompt);
    let reply = '';
    for (const word of `#${this.transcript.length} ${prompt.toUpperCase()}`.split(' ')) {
      await DG.delay(100);
      if (options?.signal?.aborted)
        throw options.signal.reason;
      reply += `${word} `;
      yield `${word} `;
    }
    return {text: reply.trim(), usage: {inputTokens: prompt.length, outputTokens: reply.length}};
  }
}

class EchoEngine extends DG.AIEngine {
  id = 'sample-echo';
  name = 'Echo (ApiSamples)';
  async available() {return false;} // never the session's default engine
  async chat(options = {}) {return new EchoChat(options);}
}

grok.ai.registerEngine(new EchoEngine());

const echo = await grok.ai.getEngine('sample-echo');
const {text} = await echo.prompt('hello world');
grok.shell.info(text);
let streamed = '';
for await (const chunk of grok.ai.stream('streamed word by word', {engine: 'sample-echo'}))
  streamed += chunk;
grok.shell.info(streamed);
