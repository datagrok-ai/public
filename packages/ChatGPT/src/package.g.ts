import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//input: string question 
export async function deepDemo(question: string) : Promise<void> {
  await PackageFunctions.deepDemo(question);
}

//tags: autostart
export function autostart() : void {
  PackageFunctions.autostart();
}

//tags: searchProvider
//output: dynamic result
export function askHelpLLMProvider() : any {
  return PackageFunctions.askHelpLLMProvider();
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

//tags: searchProvider
//output: dynamic result
export function smartChainExecutionProvider() : any {
  return PackageFunctions.smartChainExecutionProvider();
}

//input: string userGoal 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function getExecutionPlan(userGoal: string) : Promise<string> {
  return await PackageFunctions.getExecutionPlan(userGoal);
}

//input: string prompt 
//input: list<string> searchPatterns 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function fuzzyMatch(prompt: string, searchPatterns: string[]) : Promise<string> {
  return await PackageFunctions.fuzzyMatch(prompt, searchPatterns);
}

//input: string prompt 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function askDocumentationCached(prompt: string) : Promise<string> {
  return await PackageFunctions.askDocumentationCached(prompt);
}
