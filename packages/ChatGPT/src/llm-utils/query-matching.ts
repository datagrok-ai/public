/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChatGPTPromptEngine, GeminiPromptEngine, PromptEngine} from '../prompt-engine/prompt-engine';
import {_package} from '../package';
import {dartLike, fireAIAbortEvent} from '../utils';
import {AIPanelFuncs, MessageType} from './panel';
import {generateAISqlQueryWithTools} from './sql-tools';
import {ModelType} from './openAI-client';
import OpenAI from 'openai';


/// Prompt API: https://developer.chrome.com/docs/ai/prompt-api
/// Enable all Gemini Nano-related flags here:
/// chrome://flags/#prompt-api-for-gemini-nano

let _progress: DG.TaskBarProgressIndicator | null = null;

export interface QueryMatchResult {
  searchPattern: string;
  parameters: Record<string, string>;
  confidence: number;
  suggestedConnection?: string | null;
  reasoning?: string;
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

type PatternQueryMatchResult = {
  searchPattern: string;
  parameters: Record<string, string>;
  confidence: number;
  suggestedConnection?: string | null;
  reasoning?: string;
};

export async function findBestMatchingQuery(
  question: string,
): Promise<QueryMatchResult | null> {
  const schema = {
    type: 'object',
    properties: {
      searchPattern: {type: 'string'},
      parameters: {type: 'object',
        properties: {},
        required: [],
        patternProperties: {
          '.*': {},
        },
        additionalProperties: false,
      },
      confidence: {type: 'number'},
      suggestedConnection: {type: ['string', 'null']},
      reasoning: {type: 'string'},
    },
    required: ['searchPattern', 'parameters', 'confidence', 'suggestedConnection', 'reasoning'],
    additionalProperties: false,
  } as const;
  const tableQueriesSearchFunctions = DG.Func.find({meta: {searchPattern: null}, returnType: 'dataframe'})
    .filter((f) => f.options['searchPattern']);

  const searchPatterns = tableQueriesSearchFunctions.map((f) => f.options['searchPattern'] as string);
  const descriptions = tableQueriesSearchFunctions.map((f) => f.description ?? 'No description');
  const connectionIds = tableQueriesSearchFunctions.map((c) => (c instanceof DG.DataQuery) ? c.connection.id ?? 'unknown' : 'unknown');
  const idSet = new Set(connectionIds.filter((id) => id !== 'unknown'));
  const allConnections = await grok.dapi.connections.list();
  const validConnections = allConnections.filter((conn) => idSet.has(conn.id));
  const connections = connectionIds.map((id) => {
    const conn = validConnections.find((c) => c.id === id);
    return conn ? conn.nqName : 'unknown';
  });

  const systemPrompt = `
You are a precise semantic query matcher with strict validation rules. Your goal is to match user requests to database queries ONLY when you are highly confident.

CRITICAL RULES TO PREVENT HALLUCINATION:
1. NEVER fabricate or guess information not present in the user's request
2. ONLY match if the user request clearly relates to the query's purpose
3. Be conservative with confidence scores - use 0.8+ only for clear matches
4. Parameters MUST be explicitly mentioned or strongly implied in the user's request
5. Connection names are critical - mismatched connections should significantly lower confidence
6. If unsure between multiple patterns, choose the more specific one but lower the confidence

CONFIDENCE SCORING GUIDE:
- 0.9-1.0: Perfect semantic match, all parameters clear, connection aligns
- 0.7-0.89: Good match but minor ambiguity in parameters or connection
- 0.5-0.69: Weak match, significant uncertainty or missing context
- Below 0.5: Poor match, should not be used

YOUR TASK:

STEP 1 — Evaluate each query pattern:
  - Does the user's intent semantically match the query description?
  - Are the required parameters explicitly or clearly implied in the user's request?
  - Does the connection name (if not "unknown") match the user's context?
  - If connection mismatches, reduce confidence significantly

STEP 2 — Connection Inference (especially for LOW CONFIDENCE scenarios):
  - Look for ANY clues in the user's request about which data source they might want
  - This includes: explicit database/connection names, domain keywords, data types, or context
  - Be LIBERAL with connection suggestions - it's better to suggest a connection than return null
  - Even if confidence in pattern matching is low, try to infer the most likely connection
  - Match connection names flexibly (e.g., "Northwind", "northwind db", "nw" could all map to Northwind)
  - If multiple connections seem possible, choose the most likely one based on context
  - Only return null for suggestedConnection if there's truly no indication of data source

STEP 3 — Extract parameters:
  - Identify parameters in format \${parameterName} from the selected pattern
  - Extract corresponding values ONLY if explicitly stated or clearly implied
  - Use empty string if a parameter cannot be determined
  - Return as JSON object

STEP 4 — Provide reasoning:
  - In 1-2 sentences, explain why you chose this pattern and confidence level
  - Mention any concerns or ambiguities

OUTPUT REQUIREMENTS:
  - Return ONLY valid JSON matching the exact schema
  - NO markdown, NO backticks, NO code blocks, NO explanations outside JSON
  - ALL required fields must be present: searchPattern, parameters, confidence, suggestedConnection, reasoning
  - Confidence must be a number between 0 and 1
  - suggestedConnection must be a string from the connections list or null
`;

  const userPrompt = `
Available connections (nqNames):
${Array.from(new Set(connections.filter((c) => c !== 'unknown'))).join(', ')}

Available search patterns with metadata:
${searchPatterns.map((p, i) => `
Pattern #${i + 1}:
  Description: ${descriptions[i] || 'No description provided'}
  Connection: ${connections[i]}
  Pattern: ${p}
`).join('\n')}

User request:
"${question}"

Analyze the user request and return a JSON object with:
- searchPattern: the exact pattern string from above (or closest match if confidence is low)
- parameters: extracted parameter values as object
- confidence: your confidence score (0 to 1)
- suggestedConnection: infer the most likely connection from available list (be generous with inference); return null only if truly no clues exist
- reasoning: brief explanation of your decision, including why you chose/didn't choose a connection
`;

  try {
    let engine: PromptEngine;

    const isGeminiAvailable = await ('LanguageModel' in window ? LanguageModel : null)?.availability();
    if (isGeminiAvailable === 'available') {
      _package.logger.info('Using built-in Gemini model for fuzzy matching.');
      engine = GeminiPromptEngine.getInstance(geminiDownloadMonitor);
    } else {
      _package.logger.info('Using GPT engine for fuzzy matching.');
      engine = ChatGPTPromptEngine.getInstance(ModelType.Fast);
    }

    const responseText = await engine.generate(userPrompt, systemPrompt, schema);

    // Robust parsing with validation
    let result: QueryMatchResult;
    try {
      // Remove any potential markdown formatting that might slip through
      const cleanedResponse = responseText
        .replace(/^```json\s*/i, '')
        .replace(/^```\s*/i, '')
        .replace(/```\s*$/i, '')
        .trim();

      result = JSON.parse(cleanedResponse) as QueryMatchResult;

      // Validate required fields
      if (typeof result.searchPattern !== 'string')
        throw new Error('searchPattern must be a string');
      if (typeof result.parameters !== 'object' || result.parameters === null)
        throw new Error('parameters must be an object');
      if (typeof result.confidence !== 'number' || result.confidence < 0 || result.confidence > 1)
        throw new Error('confidence must be a number between 0 and 1');
      if (result.suggestedConnection !== undefined && result.suggestedConnection !== null && typeof result.suggestedConnection !== 'string')
        throw new Error('suggestedConnection must be a string or null');
      if (result.reasoning !== undefined && typeof result.reasoning !== 'string')
        throw new Error('reasoning must be a string');

      // Set defaults for optional fields if missing
      result.suggestedConnection = result.suggestedConnection ?? null;
      result.reasoning = result.reasoning ?? 'No reasoning provided';
    } catch (parseError) {
      _package.logger.error(`Failed to parse LLM response: ${parseError}. Response was: ${responseText}`);
      return null;
    }

    _package.logger.info(`Match found with confidence ${result.confidence}. Reasoning: ${result.reasoning}`);
    if (result.suggestedConnection)
      _package.logger.info(`Suggested connection: ${result.suggestedConnection}`);

    return result;
  } catch (error) {
    _package.logger.error(`Error finding best function: ${error}`);
    return null;
  }
}

export async function tableQueriesFunctionsSearchLlm(s: string): Promise<DG.Widget> {
  const tvWidgeFunc = DG.Func.find({name: 'getFuncTableViewWidget'})[0];
  if (!tvWidgeFunc)
    return DG.Widget.fromRoot(ui.divText('Power Pack plugin not installed'));
  try {
    const tableQueriesSearchFunctions = DG.Func.find({meta: {searchPattern: null}, returnType: 'dataframe'})
      .filter((f) => f.options['searchPattern']);
    const matches: PatternQueryMatchResult | null = JSON.parse(await grok.functions
      .call('ChatGPT:findMatchingPatternQuery', {prompt: s}));

    if (!matches)
      return DG.Widget.fromRoot(ui.divText('Unable to process query. Please try rephrasing.'));

    _package.logger.info(`Match result - Confidence: ${matches.confidence}, Reasoning: ${matches.reasoning || 'N/A'}`);

    if (matches.confidence > 0.7 && matches.searchPattern && matches.parameters) {
      console.log(`Function matched with high confidence ${matches.confidence}`);
      const sf = tableQueriesSearchFunctions.find((f) => {
        const pattern: string = removeTrailingQuotes(f.options['searchPattern']) ?? '';
        return pattern === matches.searchPattern;
      });

      if (sf) {
        const widget = await tvWidgeFunc.apply({func: sf, inputParams: matches.parameters});
        return widget;
      }
      _package.logger.warning(`Pattern ${matches.searchPattern} not found in available functions`);
    } else {
      console.info(`Low confidence (${matches.confidence}). ${matches.reasoning || 'No suitable match found.'}`);

      // If there's a suggested connection, provide helpful feedback
      if (matches.suggestedConnection && matches.suggestedConnection.includes(':')) {
        const nqNameSplit = matches.suggestedConnection.split(':');
        const connection = await grok.dapi.connections.filter(`namespace = "${nqNameSplit[0]}:" and shortName = "${nqNameSplit[1]}"`).first();
        if (!connection) {
          return DG.Widget.fromRoot(ui.divText(
            `Low confidence match. Your query seems related to connection "${matches.suggestedConnection}", ` +
            `but the connection was not found. Please refine your query with more specific details.`
          ));
        }
        return runQueryAIprompt(s, connection, tvWidgeFunc);
      }

      return DG.Widget.fromRoot(ui.divText(
        'No matching query or connection found. Please try rephrasing your request or check available query patterns.'
      ));
    }
  } catch (error) {
    console.error('Error calling ChatGPT:findMatchingPatternQuery', error);
    return DG.Widget.fromRoot(ui.divText(`Error during AI-powered table query search. Check console for details.`));
  }
  return DG.Widget.fromRoot(ui.divText('No matching query found via AI'));
}

function runQueryAIprompt(userPrompt: string, connection: DG.DataConnection, tvWidgetFunc: DG.Func): DG.Widget {
  const outputDiv = ui.divV([
    ui.h2(`Looks like you are trying to query the "${connection.name}" database`),
    ui.divText('Your prompt did not match any known query patterns.'),
    ui.divText('Would you like me to help you construct a query?'),
    ui.button('Generate Query', async () => {
      await runPrompt();
    })
  ], {style: {width: '100%', height: '100%'}});
  const outputWidget = DG.Widget.fromRoot(outputDiv);
  const abort = () => {
    fireAIAbortEvent();
  };


  const dummyAIPanelFuncs: AIPanelFuncs<MessageType> = {
    addUserMessage: (_aiMsg: any, _msg: string) => {},
    addAIMessage: (_aiMsg: any, _title: string, _msg: string) => {},
    addEngineMessage: (_aiMsg: any) => {},
    addUiMessage: (msg: string, fromUser: boolean, _messageOptions?: any) => {
      if (!fromUser)
        ui.setUpdateIndicator(outputDiv, true, msg.substring(0, 40) + (msg.length > 40 ? '...' : ''), abort);
    },
    addConfirmMessage: (_msg?: string) => Promise.resolve(true),
  };

  async function runPrompt() {
    try {
      ui.setUpdateIndicator(outputDiv, true, 'Generating query...', abort);
      const schemas = await grok.dapi.connections.getSchemas(connection);
      const defaultSchema = schemas.includes('public') ? 'public' : schemas.includes(connection.name) ? connection.name : schemas[0];
      const sqlResult = await generateAISqlQueryWithTools(
        userPrompt,
        connection.id!,
        defaultSchema,
        {
          aiPanel: dummyAIPanelFuncs,
          modelName: ModelType['Deep Research'],
          disableVerbose: true
        }
      );
      ui.setUpdateIndicator(outputDiv, false);
      if (sqlResult && sqlResult.trim().length > 0) {
        const queryFunc = connection.query(userPrompt, sqlResult);
        const outWidget = await tvWidgetFunc.apply({func: queryFunc, inputParams: {}});
        if (!outWidget)
          throw new Error('Failed to generate widget from AI query');
        // replace the output div content with the result widget
        outputDiv.replaceWith(outWidget.root);
        setTimeout(() => {
          const plusIcon = outWidget.root.querySelector('.d4-dialog-header .grok-icon.fal.fa-plus');
          if (!plusIcon)
            return;
          const editIcon = ui.icons.edit(() => {
            setTimeout(() => {
              const qf = connection.query(userPrompt, sqlResult);
              const v = DG.DataQueryView.create(qf);
              grok.shell.addView(v);
              // grok.shell.v = av;
            });
          }, 'Edit query');
          dartLike(editIcon.style).set('fontSize', '14px').set('width', '20px');
          plusIcon.parentNode?.insertBefore(editIcon, plusIcon);
        }, 1000);
      } else
        grok.shell.error('No SQL query was generated from the prompt. Please try rephrasing your request.');
    } catch (error: any) {
      ui.setUpdateIndicator(outputDiv, false);
      grok.shell.error(`Error during AI query generation`);
      console.error('Error during AI query generation:', error);
    }
  };
  return outputWidget;
}

function removeTrailingQuotes(s: string): string {
  let ms = s;
  if (s.startsWith('"') || s.startsWith('\''))
    ms = ms.substring(1);
  if (s.endsWith('"') || s.endsWith('\''))
    ms = ms.substring(0, ms.length - 1);
  return ms;
}
