import {LLMCredsManager} from '../llm-utils/creds';
import {OpenAIHelpClient} from '../llm-utils/openAI-client';
export interface PromptEngine {
  generate(prompt: string, system: string): Promise<string>;
}

export class ChatGPTPromptEngine implements PromptEngine {
  constructor(
    private apiKey: string,
    private model = 'gpt-4o-mini',
    private temperature = 0.0
  ) {}

  async generate(prompt: string, system: string): Promise<string> {
    const res = OpenAIHelpClient.getInstance();
    return await res.generalPromptCached(this.model, system, prompt);
  }
}

export class GeminiPromptEngine implements PromptEngine {
  constructor(
    private schema: any,
    private monitor?: CreateMonitorCallback,
  ) {}

  async generate(prompt: string, system: string): Promise<string> {
    const session = await LanguageModel.create({
      monitor: this.monitor,
      initialPrompts: [{role: 'system', content: system}],
    });

    const output = await session.prompt(prompt, {
      responseConstraint: this.schema,
    });

    return output;
  }
}
