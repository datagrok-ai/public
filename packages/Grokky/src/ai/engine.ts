import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {ClaudeModel, ClaudeRuntimeClient, FinalEvent, ImageAttachment} from '../claude/runtime-client';

export class ClaudeEngine extends DG.AIEngine {
  id = 'claude';
  name = 'Claude (Grokky)';

  private client = ClaudeRuntimeClient.getInstance();

  available(): Promise<boolean> {
    return this.client.discover();
  }

  async models(): Promise<string[]> {
    return Object.values(ClaudeModel);
  }

  async chat(options: DG.AIChatOptions = {}): Promise<DG.AIChat> {
    return new ClaudeChat(this.client, options);
  }
}

/** The runtime keeps the transcript under the session id; it has no system prompt field, so
 * `systemPrompt` goes in front of the first message. */
class ClaudeChat extends DG.AIChat {
  private sessionId = `grok-ai-${crypto.randomUUID()}`;
  private started = false;

  constructor(private client: ClaudeRuntimeClient, private options: DG.AIChatOptions) {
    super();
  }

  protected async* respond(prompt: string, options: DG.AIRunOptions = {}): DG.AIStream {
    const {client, sessionId, options: {systemPrompt, model}} = this;
    const images = await Promise.all((options.attachments ?? []).map(ClaudeChat.encode));
    await client.ensureConnected();
    const chunks: string[] = [];
    let final = null as FinalEvent | null;
    let failure: any = null;
    let wake = () => {};
    const fail = (e: any) => { failure = e; wake(); };
    const onAbort = () => fail(options.signal!.reason);
    const forSession = <T extends {sessionId: string}>(source: rxjs.Observable<T>, handler: (e: T) => void) =>
      source.subscribe((e) => { if (e.sessionId === sessionId) handler(e); });
    const subs = [
      forSession(client.onChunk, (e) => { chunks.push(e.content); wake(); }),
      forSession(client.onFinal, (e) => { final = e; wake(); }),
      forSession(client.onError, (e) => fail(new Error(e.message))),
      forSession(client.onSessionReset, () => fail(new Error('Claude runtime: the session was lost'))),
      forSession(client.onAuthRequired, () => fail(new Error('Claude runtime: authentication required'))),
      client.onClose.subscribe(() => fail(new Error('Claude runtime: connection lost'))),
    ];
    options.signal?.addEventListener('abort', onAbort);
    try {
      if (options.signal?.aborted)
        throw options.signal.reason;
      client.send(sessionId, this.started || !systemPrompt ? prompt : `${systemPrompt}\n\n${prompt}`, {
        systemPromptMode: 'none', model: model as ClaudeModel | undefined, outputSchema: options.schema, images,
      });
      while (!final || chunks.length) {
        if (chunks.length)
          yield chunks.shift()!;
        else if (failure)
          throw failure;
        else
          await new Promise<void>((resolve) => wake = resolve);
      }
      this.started = true;
      const m = final.metrics;
      return {
        text: final.content,
        structuredOutput: final.structured_output,
        usage: {
          inputTokens: m?.inputTokens ?? undefined,
          outputTokens: m?.outputTokens ?? undefined,
          costUsd: m?.costUsd ?? undefined,
        },
      };
    } finally {
      for (const s of subs)
        s.unsubscribe();
      options.signal?.removeEventListener('abort', onAbort);
      if (!final)
        client.abort(sessionId);
    }
  }

  /** The runtime forwards images only. */
  private static encode(a: DG.AIAttachment): Promise<ImageAttachment> {
    if (a.type !== 'image')
      throw new Error(`Claude engine: ${a.type} attachments are not supported`);
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve({
        mediaType: a.data.type as ImageAttachment['mediaType'], data: (reader.result as string).split(',')[1],
      });
      reader.onerror = () => reject(reader.error);
      reader.readAsDataURL(a.data);
    });
  }
}
