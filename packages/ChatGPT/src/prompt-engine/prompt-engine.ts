import {ChatModel} from 'openai/resources/shared';
import {ModelType, OpenAIClient} from '../llm-utils/openAI-client';
import {JsonSchema} from './interfaces';
export interface PromptEngine {
  generate(prompt: string, system: string, schema?: JsonSchema): Promise<string>;
}

export class ChatGPTPromptEngine implements PromptEngine {
  private static instance: ChatGPTPromptEngine | null = null;

  private constructor(
    private model = ModelType['Deep Research'],
    private temperature = 0.0
  ) {}

  public static getInstance(model?: ChatModel, temperature?: number): ChatGPTPromptEngine {
    return ChatGPTPromptEngine.instance ??= new ChatGPTPromptEngine(model, temperature);
  }

  async generate(prompt: string, system: string, schema?: JsonSchema): Promise<string> {
    const res = OpenAIClient.getInstance();
    return await res.generalPromptCached(this.model, system, prompt, schema);
  }
}

export class GeminiPromptEngine implements PromptEngine {
  private static instance: GeminiPromptEngine | null = null;
  private session: LanguageModel | null = null;

  private constructor(
    private readonly monitor?: CreateMonitorCallback,
  ) {}

  public static getInstance(monitor?: CreateMonitorCallback): GeminiPromptEngine {
    return GeminiPromptEngine.instance ??= new GeminiPromptEngine(monitor);
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
      responseConstraint: schema ?? undefined,
    });
  }
}
