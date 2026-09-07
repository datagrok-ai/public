// Streams a response. Cancel with an AbortSignal, or by leaving the loop.
// The generator's return value is the same AIResult prompt() resolves to.

const ac = new AbortController();
const output = ui.divText('');
grok.shell.newView('AI: stream', [ui.button('Stop', () => ac.abort()), output]);

const s = grok.ai.stream('Count to five.', {signal: ac.signal});
let step;
while (!(step = await s.next()).done)
  output.textContent += step.value;
grok.shell.info(`${step.value.usage.outputTokens ?? '?'} output tokens`);
