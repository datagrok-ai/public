export interface PromptEngine {
  generate(prompt: string, system?: string): Promise<string>;
}

export class ChatGPTPromptEngine implements PromptEngine {
  private url = 'https://api.openai.com/v1/chat/completions';

  constructor(
    private apiKey: string,
    private model = 'gpt-4o-mini',
    private temperature = 0.0
  ) {}

  async generate(prompt: string, system?: string): Promise<string> {
    const response = await fetch(this.url, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        Authorization: `Bearer ${this.apiKey}`,
      },
      body: JSON.stringify({
        model: this.model,
        temperature: this.temperature,
        messages: [
          { role: 'system', content: system ?? 'You are a helpful assistant.' },
          { role: 'user', content: prompt },
        ],
      }),
    });

    if (!response.ok) throw new Error(`ChatGPT API error: ${response.statusText}`);
    const result = await response.json();
    return result.choices[0].message.content;
  }
}

// Placeholder for Gemini integration
export class GeminiPromptEngine implements PromptEngine {
  async generate(prompt: string, system?: string): Promise<string> {
    if (!('ai' in window)) throw new Error('Gemini not available in this environment');
    const session = await (window as any).ai.createTextSession({ system });
    return session.prompt(prompt);
  }
}
