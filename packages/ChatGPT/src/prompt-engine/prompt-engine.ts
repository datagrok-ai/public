import {LLMCredsManager} from '../llm-utils/creds';
import {OpenAIHelpClient} from '../llm-utils/openAI-client';
export interface PromptEngine {
  generate(prompt: string, system: string, schema?: { [key: string]: unknown }): Promise<string>;
}

export class ChatGPTPromptEngine implements PromptEngine {
  constructor(
    private apiKey: string,
    private model = 'gpt-4o-mini',
    private temperature = 0.0
  ) {}

  async generate(prompt: string, system: string, schema?: { [key: string]: unknown }): Promise<string> {
    const res = OpenAIHelpClient.getInstance();
    return await res.generalPromptCached(this.model, system, prompt, schema);
  }
}

export class GeminiPromptEngine implements PromptEngine {
  private session: LanguageModel | null = null;

  constructor(
    private readonly schema: any,
    private readonly monitor?: CreateMonitorCallback,
  ) {}

  private async initSession(system: string) {
    if (!this.session) {
      this.session = await LanguageModel.create({
        monitor: this.monitor,
        initialPrompts: [{role: 'system', content: system}],
      });
    }
  }

  async generate(prompt: string, system: string): Promise<string> {
    await this.initSession(system);
    const output = await this.session!.prompt(prompt, {
      responseConstraint: this.schema,
    });
    return output;
  }
}
