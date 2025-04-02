import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type ChatGptFuncParams = { [name: string]: { type: string, description: string, default: any } };
type Message = 
  | { role: 'system' | 'user' | 'assistant'; content: string; }
  | { role: 'assistant'; content: null; function_call: null; }
  | { role: 'function'; name: string; content: string; };

export class ChatGptAssistant {
  private apiKey: string;
  private model: string;
  private temperature: number;
  private url: string = 'https://api.openai.com/v1/chat/completions';

  constructor(apiKey: string, model: string = 'gpt-4o-mini', temperature: number = 0.1) {
    this.apiKey = apiKey;
    this.model = model;
    this.temperature = temperature;
  }

  private async chatGpt(chatRequest: any): Promise<any> {
    const response = await fetch(this.url, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${this.apiKey}`,
      },
      body: JSON.stringify(chatRequest),
    });

    if (!response.ok)
      throw new Error('Failed to communicate with ChatGPT');

    return (await response.json()).choices[0];
  }

  private async executeFunction(functionName: string, parameters: any): Promise<any> {
    const func = DG.Func.find({ name: functionName })[0];
    if (!func)
      throw new Error(`Function ${functionName} not found`);
    return await func.apply(parameters);
  }

  private getType(type: string): string {
    switch (type) {
      case DG.TYPE.STRING: return 'string';
      case DG.TYPE.INT: return 'integer';
      case DG.TYPE.FLOAT: return 'number';
      case DG.TYPE.BOOL: return 'boolean';
      default: return 'object';
    }
  }

  private getProperties(f: DG.Func): ChatGptFuncParams {
    const props: ChatGptFuncParams = {};
    for (const p of f.inputs) {
      props[p.name] = {
        type: this.getType(p.propertyType),
        description: p.description ?? '',
        default: p.options['default'],
      };
    }
    return props;
  }

  private applyDefaultValues(params: any, properties: ChatGptFuncParams): any {
    const filledParams: any = {};
    for (const key in properties) {
      let value = params[key] !== undefined ? params[key] : properties[key].default;
      // If value is a string and wrapped in extra quotes, parse it safely
      if (typeof value === 'string' && value.startsWith('"') && value.endsWith('"')) {
        try {
          value = JSON.parse(value);
        } catch {

        }
      }
      filledParams[key] = value;
    }
    return filledParams;
  }

  public async askMultiStep(question: string) {
    const packages = ['Chem', 'Chembl'];
  
    const functions = packages.flatMap((pkg) =>
      DG.Func.find({ package: pkg })
        .filter(
          (f) =>
            (pkg === 'Chembl' && f.name !== 'detectChemblId') ||
            (pkg === 'Chem' && f.name === 'ChemistryGasteigerPartialCharges')
        )
        .map((f) => ({
          name: f.name,
          description: (f.nqName || f.package.name) + ": " + (f.description || f.options['description'] || f.friendlyName || f.name),
          parameters: { type: 'object', properties: this.getProperties(f) },
        }))
    );
  
    const messages: Message[] = [
      {
        role: 'system',
        content: `You are an AI assistant designed to execute JavaScript functions.
        - Your response MUST ALWAYS be structured as a function call with a name and arguments.
        - DO NOT provide free-text explanations.
        - If parameters are missing, use their default values.
        `,
      },
      { role: 'user', content: question },
    ];
  
    let intermediateResult;
    let isComplete = false;
    let executedFunction;
    let functionArguments = {};
    let executionCount = 0;
    const maxExecutions = 2;
  
    while (!isComplete && executionCount < maxExecutions) {
      const result = await this.chatGpt({
        model: this.model,
        messages: messages,
        functions: functions,
        function_call: "auto",
      });
  
      if (result.message.function_call) {
        const { name, arguments: args } = result.message.function_call;
        const parsed = JSON.parse(args);
        executedFunction = name;
  
        const func = DG.Func.find({ name })[0];
        if (!func) throw new Error(`Function "${name}" not found`);
  
        const properties = this.getProperties(func);
        const parameters = this.applyDefaultValues(parsed, properties);
        functionArguments = parameters;
  
        try {
          const functionResult = await this.executeFunction(name, parameters);
          intermediateResult = functionResult;
          const isString = typeof functionResult === 'string';
          const isHtml = functionResult instanceof HTMLElement;
  
          messages.push(
            { role: 'assistant', content: null, function_call: result.message.function_call },
            { role: 'function', name: name, content: isString ? functionResult : isHtml ? 'HTML Object generated' : "Object generated" }
          );
  
          executionCount++;
          if (executionCount >= maxExecutions) {
            isComplete = true;
          }
  
        } catch (error: any) {
          return {
            function: name,
            arguments: parameters,
            result: `Error: ${error.message}`,
          };
        }
      } else {
        isComplete = true;
      }
    }
  
    if (intermediateResult instanceof HTMLElement && intermediateResult.style.backgroundImage) {
      const backgroundImageUrl = intermediateResult.style.backgroundImage;
      let base64Data = backgroundImageUrl.split('data:image/png;base64,')[1];
      base64Data = base64Data.slice(0, -2);
      intermediateResult = ui.image(`data:image/png;base64,${base64Data}`, 200, 300);
    }
  
    return {
      function: executedFunction,
      arguments: functionArguments,
      result: intermediateResult,
    };
  }  

  public async ask(question: string): Promise<string> {
    const result = await this.chatGpt({
      model: this.model,
      messages: [{ role: 'user', content: question }],
      max_tokens: 100,
      temperature: this.temperature,
    });

    return result.message.content;
  }
}
