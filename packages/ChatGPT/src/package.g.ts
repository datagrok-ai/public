import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//tags: autostart
export function autostart() : void {
  PackageFunctions.autostart();
}

//input: string question 
//output: string result
export async function ask(question: string) : Promise<string> {
  return await PackageFunctions.ask(question);
}

//input: string question 
//output: string result
export async function askFun(question: string) : Promise<string> {
  return await PackageFunctions.askFun(question);
}
