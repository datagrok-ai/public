import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//tags: autostart
//meta.role: autostart
export function autostart() : void {
  PackageFunctions.autostart();
}

//output: dynamic result
//meta.role: searchProvider
export function combinedLLMSearchProvider() : any {
  return PackageFunctions.combinedLLMSearchProvider();
}

//name: Help
//description: Get answers from AI assistant based on Datagrok documentation and public code.
//input: string prompt 
//input: string sessionId { optional: true }
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: If the user is asking questions about how to do something, how to write the code on platform, how to execute tasks, or any other questions related to Datagrok platform functionalities and capabilities. The tone of the prompt should generally sound like "how do I do this" / "what is this". for example, "what sequence notations are supported?
export async function askHelpLLMProvider(prompt: string, sessionId?: string) : Promise<any> {
  return await PackageFunctions.askHelpLLMProvider(prompt, sessionId);
}

//description: Run the Grokky latency/accuracy benchmark suite (files/benchmark/suite.yaml) and download a JSON + Markdown report tagged with the given label. Run after logging in; open no special view.
//input: string label { description: Config label for this run, e.g. baseline / medium-effort }
//input: int reps { optional: true; description: Repetitions per prompt (default 3) }
//output: string result
export async function runBenchmark(label: string, reps?: number) : Promise<string> {
  return await PackageFunctions.runBenchmark(label, reps);
}

//description: Diff two saved benchmark runs (by label) into a Markdown delta report and download it.
//input: string labelA { description: Baseline run label }
//input: string labelB { description: Comparison run label }
//output: string result
export async function compareBenchmarks(labelA: string, labelB: string) : Promise<string> {
  return await PackageFunctions.compareBenchmarks(labelA, labelB);
}

//name: Execute
//description: Plans and executes function steps to achieve needed results
//input: string prompt 
//input: string sessionId { optional: true }
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. This relates to functions that analyse or mutate data, not get it. for example, adme properties of CHEMBL1234, enumerate some peptide, etc... Also, if the tone of the prompt sounds like "Do something to something", use this function
export async function smartChainExecutionProvider(prompt: string, sessionId?: string) : Promise<any> {
  return await PackageFunctions.smartChainExecutionProvider(prompt, sessionId);
}

//name: Query
//description: Tries to find a query which has the similar pattern as the prompt user entered and executes it
//input: string prompt 
//input: string sessionId { optional: true }
//output: widget result
//meta.role: aiSearchProvider
//meta.useWhen: if the prompt suggest that the user is looking for a data table result and the prompt resembles a query pattern. for example, "bioactivity data for shigella" or "compounds similar to aspirin" or first 100 chembl compounds. there should be some parts of user prompt that could match parameters in some query, like shigella, aspirin, first 100 etc. Always use this function when user wants to get the data without any further processing or calculating
export async function llmSearchQueryProvider(prompt: string, sessionId?: string) : Promise<any> {
  return await PackageFunctions.llmSearchQueryProvider(prompt, sessionId);
}

//input: string prompt 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function findMatchingPatternQuery(prompt: string) : Promise<string> {
  return await PackageFunctions.findMatchingPatternQuery(prompt);
}

//input: view view 
//input: string connectionID 
//input: dynamic queryEditorRoot 
//input: dynamic setAndRunFunc 
//output: bool result
export async function setupAIQueryEditor(view: DG.ViewBase, connectionID: string, queryEditorRoot: any, setAndRunFunc: any) : Promise<boolean> {
  return await PackageFunctions.setupAIQueryEditor(view, connectionID, queryEditorRoot, setAndRunFunc);
}

//description: List the catalogs available on this connection
//input: view view 
//output: string result
//meta.viewType: DataQueryView
export async function listDbCatalogs(view: any) : Promise<string> {
  return await PackageFunctions.listDbCatalogs(view);
}

//description: List schemas of a catalog (defaults to the connection default catalog)
//input: view view 
//input: string catalogName { optional: true }
//output: string result
//meta.viewType: DataQueryView
export async function listDbSchemas(view: any, catalogName?: string) : Promise<string> {
  return await PackageFunctions.listDbSchemas(view, catalogName);
}

//description: List tables of a schema with row counts
//input: view view 
//input: string schemaName 
//input: string catalogName { optional: true }
//output: string result
//meta.viewType: DataQueryView
export async function listDbTables(view: any, schemaName: string, catalogName?: string) : Promise<string> {
  return await PackageFunctions.listDbTables(view, schemaName, catalogName);
}

//description: Detailed column info (types, comments, ranges, sample values) for the given tables. Table refs: catalog.schema.table, schema.table, or table
//input: view view 
//input: string tables { description: Comma-separated table references to describe }
//output: string result
//meta.viewType: DataQueryView
export async function getDbTableDetails(view: any, tables: string) : Promise<string> {
  return await PackageFunctions.getDbTableDetails(view, tables);
}

//description: Foreign-key relationships involving the given tables — use to build correct JOINs
//input: view view 
//input: string tables { description: Comma-separated table references }
//output: string result
//meta.viewType: DataQueryView
export async function listDbJoins(view: any, tables: string) : Promise<string> {
  return await PackageFunctions.listDbJoins(view, tables);
}

//description: Test-execute a SELECT (auto-LIMITed) and report row count, columns, and a sample row. Use to validate SQL before set_query_and_run
//input: view view 
//input: string sql { description: The SQL to test }
//input: string description { description: One line describing what the query does }
//output: string result
//meta.viewType: DataQueryView
export async function getSqlTestResult(view: any, sql: string, description: string) : Promise<string> {
  return await PackageFunctions.getSqlTestResult(view, sql, description);
}

//input: string dbName { choices: ["biologics","chembl"] }
export async function moveMetaToDB(dbName: string) : Promise<void> {
  await PackageFunctions.moveMetaToDB(dbName);
}

//name: indexDatabaseSchema
export async function indexDatabaseSchema() : Promise<void> {
  await PackageFunctions.indexDatabaseSchema();
}
