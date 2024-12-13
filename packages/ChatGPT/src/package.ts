/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

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
    model: 'gpt-4', // Or whichever model you want to use
    messages: [{ role: 'user', content: question }],
    max_tokens: 100,
    temperature: temperature
  });

  return result.message.content;
}

async function executeFunction(functionName: string, parameters: any) {
  const func = DG.Func.find({name: functionName})[0];
  if (func) {
    const result = await func.apply(parameters);
    return result;
  }
  throw new Error(`Function ${functionName} not found`);
}

//input: string question
//output: object result
export async function askMultiStep(question: string) {
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
    const props: ChatGptFuncParams = {};
    for (const p of f.inputs) {
      props[p.name] = {
        type: getType(p.propertyType),
        description: p.description
      };
    }
    return props;
  }

  const packages = ['Chem', 'Chembl'];

  const functions = packages.flatMap((pkg) =>
    DG.Func.find({ package: pkg })
      .filter((f) => (pkg === 'Chembl' && f.name !== 'detectChemblId') || (pkg === 'Chem' && f.name === 'ChemistryGasteigerPartialCharges'))
      .map((f) => ({
        name: f.name,
        description: f.description || '',
        parameters: {
          type: 'object',
          properties: getProperties(f),
        },
      }))
  );

  const messages: Message[] = [
    { role: 'system', content: 'You are a helpful assistant that can call JavaScript functions when needed to calculate results.' },
    { role: 'user', content: question },
  ];

  let intermediateResult;
  let isComplete = false;
  let resultLog;

  while (!isComplete) {
    const result = await chatGpt({
      model: 'gpt-4',
      messages: messages,
      functions: functions,
      function_call: "auto",
    });

    if (result.message.function_call) {
      const { name, arguments: args } = result.message.function_call;
      const parameters = JSON.parse(args);

      try {
        const functionResult = await executeFunction(name, parameters);
        intermediateResult = functionResult;
        const isString = typeof functionResult === 'string';
        const isHtml = functionResult instanceof HTMLElement;

        messages.push(
          { role: 'assistant', content: null, function_call: result.message.function_call },
          { role: 'function', name: name, content: isString ? functionResult : isHtml ? 'HTML Object generated' : "Object generated" }
        );

      } catch (error: any) {
        console.log(`Error executing function "${name}": ${error.message}`);
        return `Error: ${error.message}`;
      }
    } else {
      isComplete = true;
      if (result.message.content)
        resultLog = result.message.content;
    }
  }

  if (intermediateResult instanceof HTMLElement && intermediateResult.style.backgroundImage) {
    const backgroundImageUrl = intermediateResult.style.backgroundImage;
    let base64Data = backgroundImageUrl.split('data:image/png;base64,')[1];
    base64Data = base64Data.slice(0, -2);
    return ui.image(`data:image/png;base64,${base64Data}`, 200, 300);
  } else if (typeof intermediateResult === 'string') {
    return ui.divText(resultLog);
  }

  return resultLog;
}