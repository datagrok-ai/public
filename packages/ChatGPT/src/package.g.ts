import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: init
//tags: init
//output: dynamic result
export async function init() {
  return PackageFunctions.init();
}

//name: chatGpt
//input: dynamic chatRequest 
//output: dynamic result
export async function chatGpt(chatRequest: any) {
  return PackageFunctions.chatGpt(chatRequest);
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
