import * as DG from 'datagrok-api/dg';

// Chrome's Prompt API global (Chrome 138+); absent in other browsers and in the Node worker
declare const LanguageModel: any;

const INPUTS = [{type: 'text'}, {type: 'image'}];

export class GemmaEngine extends DG.AIEngine {
  id = 'gemma';
  name = 'Gemma (Chrome built-in)';

  async available(): Promise<boolean> {
    return typeof LanguageModel !== 'undefined' &&
      await LanguageModel.availability({expectedInputs: INPUTS}) === 'available';
  }

  async chat(options: DG.AIChatOptions = {}): Promise<DG.AIChat> {
    if (!(await this.available()))
      throw new Error('Gemma engine: the built-in model is not downloaded; call download() first');
    return new GemmaChat(await LanguageModel.create({
      initialPrompts: options.systemPrompt ? [{role: 'system', content: options.systemPrompt}] : [],
      expectedInputs: INPUTS,
    }));
  }

  async download(): Promise<void> {
    const state = typeof LanguageModel === 'undefined' ? 'unsupported browser' :
      await LanguageModel.availability({expectedInputs: INPUTS});
    if (state === 'available')
      return;
    if (state !== 'downloadable' && state !== 'downloading')
      throw new Error(`Built-in model: ${state}`);
    const pi = DG.TaskBarProgressIndicator.create('Downloading the built-in model');
    try {
      const session = await LanguageModel.create({
        expectedInputs: INPUTS,
        monitor: (m: any) => m.addEventListener('downloadprogress',
          (e: any) => pi.update(e.loaded * 100, `${Math.round(e.loaded * 100)}%`)),
      });
      session.destroy();
    } finally {
      pi.close();
    }
  }
}

class GemmaChat extends DG.AIChat {
  constructor(private session: any) {
    super();
  }

  protected async* respond(prompt: string, options: DG.AIRunOptions = {}): DG.AIStream {
    const parts = await Promise.all((options.attachments ?? []).map(GemmaChat.encode));
    const stream = this.session.promptStreaming(
      [{role: 'user', content: [{type: 'text', value: prompt}, ...parts]}],
      {signal: options.signal, responseConstraint: options.schema});
    const reader = stream.getReader();
    let text = '';
    try {
      for (let chunk = await reader.read(); !chunk.done; chunk = await reader.read()) {
        text += chunk.value;
        yield chunk.value;
      }
    } finally {
      reader.cancel().catch(() => {});
    }
    return {text, structuredOutput: options.schema ? JSON.parse(text) : undefined, usage: {}};
  }

  close(): void {
    this.session.destroy();
  }

  private static async encode(a: DG.AIAttachment): Promise<{type: string, value: any}> {
    if (a.type === 'image')
      return {type: 'image', value: a.data};
    if (a.data.type.startsWith('text/'))
      return {type: 'text', value: `${a.title ?? 'Document'}:\n${await a.data.text()}`};
    throw new Error(`Gemma engine: ${a.data.type} documents are not supported`);
  }
}
