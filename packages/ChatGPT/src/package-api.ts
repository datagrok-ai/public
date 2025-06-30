import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('ChatGPT:Init', {});
  }

  export async function ask(question: string): Promise<any> {
    return await grok.functions.call('ChatGPT:Ask', { question });
  }

  export async function askFun(question: string): Promise<any> {
    return await grok.functions.call('ChatGPT:AskFun', { question });
  }
}
