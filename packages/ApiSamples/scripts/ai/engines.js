// The AI engines visible to grok.ai: from packages (a function with `meta.role: aiEngine`) or
// grok.ai.registerEngine(). The first available one serves calls that do not name an engine.

for (const e of await grok.ai.engines()) {
  const models = e.models ? (await e.models()).join(', ') : '';
  grok.shell.info(`${e.name} (${e.id}): ${await e.available() ? 'available' : 'not available'} ${models}`);
}
