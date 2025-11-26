/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChatGPTPromptEngine, GeminiPromptEngine, PromptEngine} from '../prompt-engine/prompt-engine';
import {modelName} from '../package';
import {_package} from '../package';
import {LLMCredsManager} from '../llm-utils/creds';


/// Prompt API: https://developer.chrome.com/docs/ai/prompt-api
/// Enable all Gemini Nano-related flags here:
/// chrome://flags/#prompt-api-for-gemini-nano

let _progress: DG.TaskBarProgressIndicator | null = null;

export interface QueryMatchResult {
  searchPattern: string;
  parameters: Record<string, string>;
  confidence: number;
}

export function geminiDownloadMonitor(m: CreateMonitor) {
  m.addEventListener('downloadprogress', (e) => {
    if (!_progress) {
      grok.shell.info('Gemini model download started. This might take a few minutes.');
      _progress = DG.TaskBarProgressIndicator.create('Downloading Gemini...');
    }
    _progress.update(e.loaded * 100, 'Downloading Gemini...');
  });
}

export async function findBestFunction(
  question: string,
  searchPatterns: string[],
  descriptions: string[],
): Promise<QueryMatchResult | null> {
  const schema = {
    searchPattern: 'string',
    parameters: 'object',
    confidence: 'number',
  };

  const connections = searchPatterns.map((pattern) => {
    const func = DG.Func.find({meta: {searchPattern: pattern}})[0];
    if (!func || !(func instanceof DG.DataQuery) || !func.query)
      return 'unknown';
    const connLine = func.query.split('\n').find((line) => line.startsWith('--connection:'));
    return connLine ? connLine.replace('--connection:', '').trim() : 'unknown';
  });

  const systemPrompt = `
You are an intelligent semantic query matcher. Your job is to find which query search pattern best matches a user's request.

IMPORTANT: pay careful attention to the keywords in description and search pattern and also the connection. DO NOT ARTIFICIALLY INFLATE CONFIDENCE.

Follow these steps carefully:
STEP 1 — Identify the most relevant structured query:
  - Find the query whose intent best matches the user's request semantically, even if the phrasing or order differs.
  - Prioritize matches if the name of the connection (except unknown) aligns with user's prompt, otherwise lower the confidence!
  - Compare using meaning, not just keywords.
  - If several are close, choose the one with the most specific match.

STEP 2 — Extract parameter values:
  - For the selected query, identify all parameters written as \${parameter}.
  - Extract the corresponding values from the user input or infer them if clearly implied.
  - Return them in JSON form.

Important instructions:
  - **Output only valid JSON.**
  - **Do not include markdown, backticks, or any extra characters.**
  - **Do not provide explanations or commentary.**
  - Your response must exactly match this structure:

  ${JSON.stringify(schema, null, 2)}
`;

  const userPrompt = `
Available search patterns:
${searchPatterns.map((p, i) => `${i + 1}. \n description: ${descriptions[i] || 'Empty'} \n Connection: ${connections[i]} \n pattern: ${p} \n`).join('\n')}

User request:
"${question}"
`;

  try {
    let engine: PromptEngine;

    const isGeminiAvailable = await ('LanguageModel' in window ? LanguageModel : null)?.availability();
    if (isGeminiAvailable === 'available') {
      _package.logger.info('Using built-in Gemini model for fuzzy matching.');
      engine = GeminiPromptEngine.getInstance(schema, geminiDownloadMonitor);
    } else {
      _package.logger.info('Using GPT engine for fuzzy matching.');
      engine = ChatGPTPromptEngine.getInstance(LLMCredsManager.getApiKey(), modelName);
    }

    const responseText = await engine.generate(userPrompt, systemPrompt);
    return JSON.parse(responseText) as QueryMatchResult;
  } catch (error) {
    _package.logger.error(`Error finding best function: ${error}`);
    return null;
  }
}

export async function tableQueriesFunctionsSearchLlm(s: string): Promise<DG.Widget> {
  const tvWidgeFunc = DG.Func.find({name: 'getFuncTableViewWidget'})[0];
  if (!tvWidgeFunc)
    return DG.Widget.fromRoot(ui.divText('Power Pack plugin not installed'));
  const tableQueriesSearchFunctions = DG.Func.find({meta: {searchPattern: null}, returnType: 'dataframe'})
    .filter((f) => f.options['searchPattern']);

  const searchPatterns = tableQueriesSearchFunctions.map((f) => f.options['searchPattern']);
  const descriptions = tableQueriesSearchFunctions.map((f) => f.description ?? 'No description');
  try {
    const matches = JSON.parse(await grok.functions
      .call('ChatGPT:fuzzyMatch', {prompt: s, searchPatterns, descriptions}));
    if (matches && matches.searchPattern && matches.parameters) {
      const sf = tableQueriesSearchFunctions.find((f) => {
        const pattern: string = removeTrailingQuotes(f.options['searchPattern']) ?? '';
        return pattern === matches.searchPattern;
      });

      if (sf) {
        const widget = await tvWidgeFunc.apply({func: sf, inputParams: matches.parameters});
        return widget;
      }
    }
  } catch (error) {
    console.error('Error calling ChatGPT:fuzzyMatch', error);
    return DG.Widget.fromRoot(ui.divText(`Error during AI-powered table query search. Check console for details.`));
  }
  return DG.Widget.fromRoot(ui.divText('No matching query found via AI'));
}

function removeTrailingQuotes(s: string): string {
  let ms = s;
  if (s.startsWith('"') || s.startsWith('\''))
    ms = ms.substring(1);
  if (s.endsWith('"') || s.endsWith('\''))
    ms = ms.substring(0, ms.length - 1);
  return ms;
}
