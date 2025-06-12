import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('Chatgpt:Init', {});
  }

  export async function ask(question: string): Promise<any> {
    return await grok.functions.call('Chatgpt:Ask', { question });
  }

  export async function askFun(question: string): Promise<any> {
    return await grok.functions.call('Chatgpt:AskFun', { question });
  }
}
