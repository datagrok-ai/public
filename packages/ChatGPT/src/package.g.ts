import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() {
  return PackageFunctions.init();
}

//name: autostart
//tags: autostart
export function autostart() {
  return PackageFunctions.autostart();
}

//name: ask
//input: string question 
//output: string result
export async function ask(question: string) {
  return PackageFunctions.ask(question);
}

//name: askFun
//input: string question 
//output: string result
export async function askFun(question: string) {
  return PackageFunctions.askFun(question);
}
