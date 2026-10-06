// A multi-turn conversation: the engine keeps the transcript. Close it when done.

const chat = await grok.ai.chat({systemPrompt: 'You are a bioassay analyst. Answer in one sentence.'});
try {
  const first = await chat.prompt('Which dose-response curve shape indicates a problem?');
  grok.shell.info(first.text);
  const second = await chat.prompt('Why that one in particular?');
  grok.shell.info(second.text);
} finally {
  chat.close();
}
