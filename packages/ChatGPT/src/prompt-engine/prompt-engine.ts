import {OpenAIHelpClient} from '../llm-utils/openAI-client';
import {JsonSchema} from './interfaces';
export interface PromptEngine {
  generate(prompt: string, system: string, schema?: JsonSchema): Promise<string>;
}

export class ChatGPTPromptEngine implements PromptEngine {
  private static instance: ChatGPTPromptEngine | null = null;

  private constructor(
    private apiKey: string,
    private model = 'gpt-4o-mini',
    private temperature = 0.0
  ) {}

  public static getInstance(apiKey: string, model?: string, temperature?: number): ChatGPTPromptEngine {
    return ChatGPTPromptEngine.instance ??= new ChatGPTPromptEngine(apiKey, model, temperature);
  }

  async generate(prompt: string, system: string, schema?: JsonSchema): Promise<string> {
    const res = OpenAIHelpClient.getInstance();
    return await res.generalPromptCached(this.model, system, prompt, schema);
  }
}

export class GeminiPromptEngine implements PromptEngine {
  private static instance: GeminiPromptEngine | null = null;
  private session: LanguageModel | null = null;

  private constructor(
    private readonly schema: any,
    private readonly monitor?: CreateMonitorCallback,
  ) {}

  public static getInstance(schema: any, monitor?: CreateMonitorCallback): GeminiPromptEngine {
    return GeminiPromptEngine.instance ??= new GeminiPromptEngine(schema, monitor);
  }

  private async initSession(system: string) {
    if (!this.session) {
      this.session = await LanguageModel.create({
        monitor: this.monitor,
        initialPrompts: [{role: 'system', content: system}],
      });
    }
  }

  async generate(prompt: string, system: string, schema?: JsonSchema): Promise<string> {
    await this.initSession(system);
    return this.session!.prompt(prompt, {
      responseConstraint: schema ?? this.schema,
    });
  }
}
