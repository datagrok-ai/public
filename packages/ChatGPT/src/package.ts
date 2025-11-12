/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {powerSearchQueryTable, processPowerSearchTableView} from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';
import { ChatGptAssistant } from './prompt-engine/chatgpt-assistant';
import { ChatGPTPromptEngine } from './prompt-engine/prompt-engine';
import { AssistantRenderer } from './prompt-engine/rendering-tools';

export const _package = new DG.Package();

let apiKey: string = '';
let model: string = 'gpt-4';
let temperature = 0.1;
const url = 'https://api.openai.com/v1/chat/completions';

//tags: init
export async function init() {
  // @ts-ignore
  apiKey = (await _package.getSettings())['apiKey'];
}

export async function chatGpt(chatRequest: any): Promise<any> {
  const response = await grok.dapi.fetchProxy(url, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${apiKey}`
    },
    body: JSON.stringify(chatRequest)
  });

  if (!response.ok)
    throw new Error('Failed to communicate with ChatGPT');

  const data = await response.json();
  return data.choices[0];
}

//input: string question
//output: string answer
export async function ask(question: string): Promise<string> {
  let result = await chatGpt({
    model: 'gpt-4o-mini', // Or whichever model you want to use
    messages: [{ role: 'user', content: question }],
    max_tokens: 100,
    temperature: temperature
  });

  return result.message.content;
}

//tags: search
//input: string question
//output: widget w
export async function askMultiStep(question: string): Promise<DG.Widget> {
  return DG.Widget.fromRoot(
    (() => {
      const container = ui.divV([]);
      
      const button = ui.button('Ask AI', async () => {
        ui.empty(container);
        const loader = ui.loader();
        container.appendChild(loader);

        try {
          const gptEngine = new ChatGPTPromptEngine(apiKey, 'gpt-4.1-mini-2025-04-14');
          const gptAssistant = new ChatGptAssistant(gptEngine);
          const { plan, result } = await gptAssistant.plan(question);
          ui.empty(container);
          container.append(AssistantRenderer.renderFullOutput(plan, result).root);
        } catch (error) {
          ui.empty(container);
          container.appendChild(ui.divText('Error while processing request.', { style: { textAlign: 'center' } }));
        }
      });
      container.style.paddingTop = '10px';
      const wrapper = ui.divV([button, container]);
      return wrapper;
    })()
  );
}