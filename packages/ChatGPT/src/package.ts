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
  // this is cheating but hey :D
  const styleElement = document.createElement('style');
  styleElement.innerHTML = `
  .power-pack-widget-host:has(.ask-ai-widget-container) {
    height: unset!important;
  }
  `;
  document.head.appendChild(styleElement);

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

async function createTableQueryWidget(func: DG.Func, inputParams: any): Promise<DG.Widget> {

  const container = ui.divV([]);
  container.style.minHeight = '30px';
  const fc = func.prepare(inputParams);
  let firstTime = true;
  Array.from(fc.inputParams.values()).forEach((i) => {
  DG.debounce(i.onChanged, 1000).subscribe(() => {
      perform();
    });
  })
  const editor = await fc.getEditor();
  editor.style.width = '100%';
  const tableContainer = ui.div([], {style: {width: '100%', height: '100%'}});
  if (editor) container.appendChild(editor);
  container.appendChild(tableContainer);
      
  const perform = async () => {
    ui.empty(tableContainer);
    ui.setUpdateIndicator(tableContainer, true);
    try {
      const resFuncCall = await fc.call();
      const views = resFuncCall.getResultViews();
      const tv = views[0] as DG.TableView;
      
      if (tv) {
        tableContainer.appendChild(tv.root);
        container.style.height = 'calc(100vh - 300px)';
        container.style.width = '100%';
        tv._onAdded();
        if (firstTime) {
          setTimeout(() => {
            processPowerSearchTableView(tv);
          }, 100);
          firstTime = false;
        }
      }
    } catch (e) {
      console.error('Error executing query function:', e);
      ui.empty(tableContainer);
      tableContainer.appendChild(ui.divText('Error executing query function'));
    }
    ui.setUpdateIndicator(tableContainer, false);
  }

  perform();

  return DG.Widget.fromRoot(container);
  
}

//tags: search
//input: string question
//output: widget w
export async function askMultiStep(question: string): Promise<DG.Widget> {
  return DG.Widget.fromRoot(
    (() => {
      const container = ui.divV([], {classes: 'ask-ai-widget-container', style: {alignItems: 'center', minHeight: '30px'}});
      
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
            container.appendChild(
              ui.divText(`Function "${response.function}" not found.`, { style: { textAlign: 'center' } })
            );
            return;
          }

          if (func.type === 'data-query' && !(func.outputs.length == 1 && func.outputs[0].propertyType === 'string')) {
            const inputParams = response.arguments;
            container.appendChild(
              (await createTableQueryWidget(func, inputParams)).root
            );
          } else {
            const result = response.result;
            if (typeof result === 'string') {
              if (grok.chem.checkSmiles(result)) {
                container.appendChild(grok.chem.drawMolecule(result, 300,300));
              } else {
                container.appendChild(ui.divText(result, { style: { textAlign: 'center' } }));
              }
            } else if (result instanceof HTMLElement) {
              container.appendChild(result);
            } else {
              container.appendChild(ui.divText('Unsupported result type: ' + typeof result , { style: { textAlign: 'center' } }));
            }
          }
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