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
//meta.useWhen: If the user is asking questions about how to do something, how to write the code on platform, how to execute tasks, or any other questions related to Datagrok platform functionalities and capabilities. The tone of the prompt should generally sound like "how do I do this" / "what is this". for example, "what sequence notations are supported?
export async function askHelpLLMProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.askHelpLLMProvider(prompt);
}

//name: Execute
//description: Plans and executes function steps to achieve needed results
//input: string prompt 
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. for example, adme properties of CHEMBL1234, enumerate some peptide, etc.. . Also, if the tone of the prompt sounds like "Do something", use this function
export async function smartChainExecutionProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.smartChainExecutionProvider(prompt);
}

//name: Query
//description: Tries to find a query which has the similar pattern as the prompt user entered and executes it
//input: string prompt 
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: if the prompt suggest that the user is looking for a data table result and the prompt resembles a query pattern. for example, "bioactivity data for shigella" or "compounds similar to aspirin" or first 100 chembl compounds. there should be some parts of user prompt that could match parameters in some query, like shigella, aspirin, first 100 etc.
export async function llmSearchQueryProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.llmSearchQueryProvider(prompt);
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

//input: string question 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function ask(question: string) : Promise<string> {
  return await PackageFunctions.ask(question);
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
//input: list<string> descriptions 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function fuzzyMatch(prompt: string, searchPatterns: string[], descriptions: string[]) : Promise<string> {
  return await PackageFunctions.fuzzyMatch(prompt, searchPatterns, descriptions);
}

//input: string prompt 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function askDocumentationCached(prompt: string) : Promise<string> {
  return await PackageFunctions.askDocumentationCached(prompt);
}

//input: string prompt 
//input: string connectionID 
//input: string schemaName 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function generateSqlQuery(prompt: string, connectionID: string, schemaName: string) : Promise<string> {
  return await PackageFunctions.generateSqlQuery(prompt, connectionID, schemaName);
}
