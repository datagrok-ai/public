/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ChatGptAssistant } from './gpt-base';
import {powerSearchQueryTable, processPowerSearchTableView} from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';

export const _package = new DG.Package();

type ChatGptFuncParams = { [name: string]: { type: string, description: string } };
type Message = 
  | { role: 'system' | 'user' | 'assistant'; content: string; }
  | { role: 'assistant'; content: null; function_call: { name: string; arguments: string; }; }
  | { role: 'function'; name: string; content: string; };


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
  const response = await fetch(url, {
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

function createTableQueryWidget(func: DG.Func, inputParams: any): DG.Widget {
  return DG.Widget.fromRoot(ui.wait(async () => {
    try {
      const fc = func.prepare(inputParams);
      const resFuncCall = await fc.call();
      const views = resFuncCall.getResultViews();
      const tv = views[0] as DG.TableView;
      const container = ui.divV([]);

      const editor = await resFuncCall.getEditor();
      if (editor) container.appendChild(editor);

      setTimeout(() => {
        processPowerSearchTableView(tv);
        tv._onAdded();
      }, 200);

      container.appendChild(tv.root);
      return container;
    } catch (e) {
      console.error('Error executing query function:', e);
      return ui.divText('Error executing query function');
    }
  }));
}

//tags: search
//input: string question
//output: widget w
export async function askMultiStep(question: string): Promise<DG.Widget> {
  return DG.Widget.fromRoot((() => {
    const container = ui.divV([]);
    const button = ui.button('Ask AI', async () => {
      ui.empty(container);
      const loader = ui.loader();
      container.appendChild(loader);

      try {
        const gptAssistant = new ChatGptAssistant(apiKey);
        const response = await gptAssistant.askMultiStep(question);

        ui.empty(container);

        const func = DG.Func.find({ name: response.function })[0];

        if (!func) {
          container.appendChild(ui.divText(`Function "${response.function}" not found.`));
          return;
        }

        if (func.type === 'data-query') {
          const inputParams = response.arguments;
          container.appendChild(createTableQueryWidget(func, inputParams).root);
        } else if (func.type === 'script') {
          container.appendChild(response.result);
        }
      } catch (error) {
        ui.empty(container);
        container.appendChild(ui.divText('Error while processing request.'));
      }
    });
    
    const wrapper = ui.divV([button, container], { style: { gap: '10px' } });

    return wrapper;
  })());
}
