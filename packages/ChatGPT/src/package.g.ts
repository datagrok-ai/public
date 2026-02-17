import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//meta.role: autostart
export function autostart() : void {
  PackageFunctions.autostart();
}

//name: setupClaudeRuntimeForTableView
export async function setupClaudeRuntimeForTableView() : Promise<void> {
  await PackageFunctions.setupClaudeRuntimeForTableView();
}

//output: dynamic result
//meta.role: searchProvider
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
//meta.useWhen: If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. This relates to functions that analyse or mutate data, not get it. for example, adme properties of CHEMBL1234, enumerate some peptide, etc... Also, if the tone of the prompt sounds like "Do something to something", use this function
export async function smartChainExecutionProvider(prompt: string) : Promise<any> {
  return await PackageFunctions.smartChainExecutionProvider(prompt);
}

//name: Query
//description: Tries to find a query which has the similar pattern as the prompt user entered and executes it
//input: string prompt 
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: if the prompt suggest that the user is looking for a data table result and the prompt resembles a query pattern. for example, "bioactivity data for shigella" or "compounds similar to aspirin" or first 100 chembl compounds. there should be some parts of user prompt that could match parameters in some query, like shigella, aspirin, first 100 etc. Always use this function when user wants to get the data without any further processing or calculating
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
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function findMatchingPatternQuery(prompt: string) : Promise<string> {
  return await PackageFunctions.findMatchingPatternQuery(prompt);
}

//input: string prompt 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function askDocumentationCached(prompt: string) : Promise<string> {
  return await PackageFunctions.askDocumentationCached(prompt);
}

//input: view view 
//input: string connectionID 
//input: dynamic queryEditorRoot 
//input: dynamic setAndRunFunc 
//output: bool result
export async function setupAIQueryEditor(view: DG.ViewBase, connectionID: string, queryEditorRoot: any, setAndRunFunc: any) : Promise<boolean> {
  return await PackageFunctions.setupAIQueryEditor(view, connectionID, queryEditorRoot, setAndRunFunc);
}

//input: string dbName { choices: ["biologics","chembl"] }
export async function moveMetaToDB(dbName: string) : Promise<void> {
  await PackageFunctions.moveMetaToDB(dbName);
}

//name: setupVectorStore
export async function setupVectorStore() : Promise<void> {
  await PackageFunctions.setupVectorStore();
}

//name: searchForSomething
export async function searchForSomething() : Promise<void> {
  await PackageFunctions.searchForSomething();
}

//name: indexDatabaseSchema
export async function indexDatabaseSchema() : Promise<void> {
  await PackageFunctions.indexDatabaseSchema();
}
