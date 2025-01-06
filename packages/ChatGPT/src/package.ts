/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

type ChatGptFuncParams = { [name: string]: { type: string, description: string } };

let apiKey: string = '';
let model: string = 'gpt-4';
let temperature = 0.7;
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
    model: 'gpt-4', // Or whichever model you want to use
    messages: [{ role: 'user', content: question }],
    max_tokens: 100,
    temperature: temperature
  });

  return result.message.content;
}

//input: string question
//output: string answer
export async function askFun(question: string): Promise<string> {

  function getType(type: string): string {
    switch (type) {
      case DG.TYPE.STRING: return 'string';
      case DG.TYPE.INT: return 'integer';
      case DG.TYPE.FLOAT: return 'number';
      case DG.TYPE.BOOL: return 'boolean';
      default: return 'object';
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
    model: 'gpt-4', // Or whichever model you want to use
    messages: [
      { role: 'system', content: 'You are a helpful assistant that can call JavaScript functions when needed.' },
      { role: 'user', content: question }
    ],
    functions: functions,
    function_call: "auto"
  });

  if (result.message.function_call) {
    const functionName = result.message.function_call.name;
    const parameters = result.message.function_call.arguments;
    const functionResult = await executeFunction(functionName, JSON.parse(parameters));
    return functionResult;
  }
  return JSON.stringify(result);
}

async function executeFunction(functionName: string, parameters: any) {
  const func = DG.Func.find({name: functionName})[0];
  if (func) {
    const result = await func.apply(parameters);
    return result;
  }
  throw new Error(`Function ${functionName} not found`);
}
