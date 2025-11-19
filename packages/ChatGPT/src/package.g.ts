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

//tags: searchProvider
//output: dynamic result
export function combinedLLMSearchProvider() : any {
  return PackageFunctions.combinedLLMSearchProvider();
}

//name: Help
//description: Get answers from DeepGROK AI assistant based on Datagrok documentation and public code.
//input: string prompt 
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: If the user is asking questions about how to do something, how to write the code on platform, how to execute tasks, or any other questions related to Datagrok platform functionalities and capabilities. for example, "what sequence notations are supported?
export async function askHelpLLMProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.askHelpLLMProvider(prompt);
}

//name: Execute
//description: Plans and executes function steps to achieve needed results
//input: string prompt 
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. for example, adme properties of CHEMBL1234
export async function smartChainExecutionProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.smartChainExecutionProvider(prompt);
}

//input: string model 
//input: string systemPrompt 
//input: string prompt 
//input: dynamic schema 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function askAIGeneralCached(model: string, systemPrompt: string, prompt: string, schema?: any) : Promise<string> {
  return await PackageFunctions.askAIGeneralCached(model, systemPrompt, prompt, schema);
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
