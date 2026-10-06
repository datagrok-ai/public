// Every turn, from any engine and any caller, is observable on grok.ai.

const subs = [
  grok.ai.onPrompt.subscribe((e) => grok.shell.info(`prompt: ${e.prompt}`)),
  grok.ai.onResult.subscribe((e) => grok.shell.info(`result: ${e.result.usage.outputTokens ?? '?'} output tokens`)),
  grok.ai.onCancel.subscribe((e) => grok.shell.info('cancelled')),
  grok.ai.onError.subscribe((e) => grok.shell.warning(e.error.message)),
];
try {
  await grok.ai.prompt('Name three colours.');
} finally {
  subs.forEach((s) => s.unsubscribe());
}
