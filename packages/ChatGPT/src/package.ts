/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ChatGptAssistant } from './prompt-engine/chatgpt-assistant';
import { ChatGPTPromptEngine } from './prompt-engine/prompt-engine';
import { AssistantRenderer } from './prompt-engine/rendering-tools';
import {getAiPanelVisibility, initAiPanel, setAiPanelVisibility} from './ai-panel';

export * from './package.g';
export const _package = new DG.Package();

type ChatGptFuncParams = { [name: string]: { type: string, description: string } };

type IChatGptResponse = {
  finish_reason?: string;
  index?: number;
  logprobs?: any;
  message?: {
    function_call?: any;
    role?: string;
    content?: string;
    refusal?: any;
  }
}

let apiKey: string = '';
let model: string = 'gpt-4';
let temperature = 0.1;
const url = 'https://api.openai.com/v1/chat/completions';


export async function chatGpt(chatRequest: any): Promise<IChatGptResponse> {
  const response = await fetch(url, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${apiKey}`
    },
    body: JSON.stringify(chatRequest)
  });

  if (!response.ok) {
    grok.shell.error('ChatGPT error: ' + response.statusText);
    throw new Error('Failed to communicate with ChatGPT');
  }

  const data = await response.json();
  return data.choices[0];
}

export async function askImpl(question: string): Promise<IChatGptResponse> {
  const request: any = {
    model: model,
    messages: [{ role: 'user', content: question }],
    max_tokens: 100,
    temperature: temperature
  }

  return await chatGpt(request);
}

//tags: search
//input: string question
//output: widget w
export async function askMultiStep(question: string): Promise<DG.Widget> {
  return DG.Widget.fromRoot(
    (() => {
      const planContainer = ui.divV([]);
      const resultContainer = ui.divV([]);
      const gptEngine = new ChatGPTPromptEngine(apiKey, 'gpt-4.1-mini-2025-04-14');
      const gptAssistant = new ChatGptAssistant(gptEngine);

      const button = ui.button('Ask AI', () => {
        ui.empty(planContainer);
        ui.empty(resultContainer);

        const planWait = ui.wait(async () => {
          const plan = await gptAssistant.plan(question);
          return AssistantRenderer.renderPlan(plan);
        });
        planContainer.appendChild(planWait);

        const resultWait = ui.wait(async () => {
          const plan = await gptAssistant.plan(question);
          const result = await gptAssistant.execute(plan);

          return AssistantRenderer.renderResult(result);
        });
        resultContainer.appendChild(resultWait);
      });

      const wrapper = ui.divV([button, planContainer, resultContainer], 'chatgpt-ask-ai-result');
      return wrapper;
    })()
  );
}

async function executeFunction(functionName: string, parameters: any) {
  const func = DG.Func.find({name: functionName})[0];
  if (func) {
    return await func.apply(parameters);
  }
  throw new Error(`Function ${functionName} not found`);
}

export class PackageFunctions {

  @grok.decorators.init()
  static async init() {
    apiKey = _package.settings['apiKey'];
  }


  @grok.decorators.autostart()
  static autostart() {
    // grok.shell.info('started')
    //
    grok.events.onViewAdded.subscribe((view) => {
      if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
        const tableView = view as DG.TableView;
        const iconFse = ui.iconSvg('ai.svg', () => setAiPanelVisibility(true), 'Ask AI');
        tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);
      }
    });

    initAiPanel();
    // Add keyboard shortcut for toggling AI panel
    document.addEventListener('keydown', (event) => {
      // Check for Ctrl+I (Ctrl key + I key)
      if (event.ctrlKey && event.key === 'i') {
        event.preventDefault(); // Prevent default browser behavior
        const isVisible = getAiPanelVisibility();
        setAiPanelVisibility(!isVisible);
      }
    });
  }


  @grok.decorators.func()
  static async ask(question: string): Promise<string> {
    let result = await askImpl(question);
    return result.message!.content!;
  }


  @grok.decorators.func()
  static async askFun(question: string): Promise<string> {
    function getType(type: string): string {
      switch (type) {
        case DG.TYPE.STRING:
          return 'string';
        case DG.TYPE.INT:
          return 'integer';
        case DG.TYPE.FLOAT:
          return 'number';
        case DG.TYPE.BOOL:
          return 'boolean';
        default:
          return 'object';
      }
    }

    function getProperties(f: DG.Func): ChatGptFuncParams {
      let props: ChatGptFuncParams = {};
      for (const p of f.inputs) {
        props[p.name] = {
          type: getType(p.propertyType),
          description: p.description
        }
      }
      return props;
    }

    const functions = DG.Func.find({package: 'Admetica'}).map((f) => {
      return {
        name: f.name,
        description: f.description,
        parameters: {
          type: "object",
          properties: getProperties(f)
        }
      }
    });

    const result = await chatGpt({
      model: 'gpt-4',
      messages: [
        {role: 'system', content: 'You are a helpful assistant that can call JavaScript functions when needed.'},
        {role: 'user', content: question}
      ],
      functions: functions,
      function_call: "auto"
    });

    if (result.message?.function_call) {
      const functionName = result.message.function_call.name;
      const parameters = result.message.function_call.arguments;
      return await executeFunction(functionName, JSON.parse(parameters));
    }
    return JSON.stringify(result);
  }
}
