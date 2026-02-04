/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AIPanelFuncs, MessageType} from './panel';
import {getAIAbortSubscription} from '../utils';
import {ModelType, LLMClient} from './LLM-client';
import {LLMCredsManager} from './creds';
import {ComparisonFilter, CompoundFilter} from 'openai/resources/shared';
import {LanguageModelV3FunctionTool, LanguageModelV3Message} from '@ai-sdk/provider';
import {findLast} from '../utils';


class ClarificationNeededError extends Error {
  constructor() {
    super('CLARIFICATION_NEEDED');
  }
}

class NoScriptGeneratedError extends Error {
  constructor() {
    super('NO_SCRIPT_GENERATED');
  }
}

type JsonPrimitive = string | number | boolean | null;
type JsonValue = JsonPrimitive | JsonObject | JsonValue[];
type JsonObject = { [key: string]: JsonValue };

function isJsonObject(value: JsonValue | object | null | undefined): value is JsonObject {
  return !!value && typeof value === 'object' && !Array.isArray(value);
}

function getStringProp(obj: JsonObject | null | undefined, key: string): string | null {
  const v = obj?.[key];
  return typeof v === 'string' ? v : null;
}

function getNumberProp(obj: JsonObject | null | undefined, key: string): number | null {
  const v = obj?.[key];
  return typeof v === 'number' ? v : null;
}

function parseJsonObject(text: string): JsonObject | null {
  try {
    const parsed = JSON.parse(text) as JsonValue;
    return isJsonObject(parsed) ? parsed : null;
  } catch {
    return null;
  }
}

function isSafeJsMemberExpression(expr: string): boolean {
  // Allow only dotted identifiers like: grok.chem, DG, ui, grok.shell.tv
  // Disallow calls, indexing, quotes, spaces, operators, etc.
  return /^[A-Za-z_$][\w$]*(\.[A-Za-z_$][\w$]*)*$/.test(expr);
}


/**
 * Generates a Datagrok script using function calling approach with Responses API
 * @throws Error if script generation fails
 * @param prompt - User's natural language description of the script
 * @param language - Programming language for the script
 * @param options - Additional options including message history and AI panel functions
 */
export async function generateDatagrokScript(
  prompt: string,
  language: DG.ScriptingLanguage,
  options: {
    oldMessages?: MessageType[],
    aiPanel?: AIPanelFuncs<MessageType>,
    disableVerbose?: boolean,
  } = {}
): Promise<string> {
  let aborted = false;

  options.aiPanel?.addUiMessage(prompt, true);

  // Create tool execution context
  const vectorStoreId = LLMCredsManager.getVectorStoreId();
  const context = new ScriptGenerationContext(vectorStoreId);
  const abortSub = getAIAbortSubscription().subscribe(() => {
    aborted = true;
    console.log('Aborting script generation as per user request');
    try {
      context.dispose();
    } catch {
      // ignore
    }
    abortSub.unsubscribe();
  });

  try {
    // Initialize OpenAI client
    // const openai = OpenAIClient.getInstance().openai;
    const langTool = LLMClient.getInstance();
    const client = langTool.aiModels.Coding;

    // System message with strict guidelines
    const systemMessage = getSystemPrompt(language, vectorStoreId);

    // Create a running input list we will add to over time (per Responses API docs)
    const input: MessageType[] = options.oldMessages ? [...options.oldMessages] : [];
    if (!input.length) {
      const systemMsg = langTool.createSystemMessage(systemMessage);
      input.push(systemMsg);
      options.aiPanel?.addEngineMessage(systemMsg);
    }
    const userMsg = langTool.createUserMessage(prompt);
    options.aiPanel?.addEngineMessage(userMsg);
    input.push(userMsg);
    // options.aiPanel?.addEngineMessage(systemMessage);
    // Define available tools
    const tools = getTools(language, vectorStoreId);

    // Function calling loop using Responses API
    let iterations = 0;
    let maxIterations = 20; // Prevent infinite loops


    while (iterations < maxIterations) {
      if (aborted) {
        grok.shell.info('SQL generation aborted by user.');
        return '';
      }
      iterations++;
      console.log(`\n=== Iteration ${iterations} ===`);

      const response = await client.doGenerate({
        prompt: input,
        tools,
        providerOptions: {
          openai: {
            ...(ModelType.Coding.startsWith('gpt-5') ? {reasoning: {effort: 'medium'}} : {}),
          }
        }
      });
      const fixedContent = response.content.map((item) => {
        if (item.type === 'tool-call') {
          return {
            ...item,
            input: typeof item.input === 'string' ? parseJsonObject(item.input) ?? item.input : item.input
          };
        }
        return item;
      });

      const formattedOutput: LanguageModelV3Message = {
        role: 'assistant',
        // @ts-ignore
        content: fixedContent
      };
      //const outputs: MessageType[] = response.content.map((item) => ({role: 'assistant', content: [item]} as MessageType));
      input.push(formattedOutput);
      options.aiPanel?.addEngineMessage(formattedOutput);
      // outputs.forEach((item) =>
      //   options.aiPanel?.addEngineMessage(item));

      // Execute function calls (if any)
      let hadToolCalls = false;
      for (const call of response.content.filter((c) => c.type === 'tool-call')) {
        hadToolCalls = true;

        const functionName = call.toolName;
        const callId = call.toolCallId;
        const rawArgs = call.input ?? '{}';
        const args = parseJsonObject(rawArgs);

        if (!functionName || !callId) {
          const errorText = 'Error: malformed function call (missing name/call_id).';
          const outputItem = langTool.createToolOutputMessage(callId ?? '', functionName ?? '', errorText);
          if (callId) {
            input.push(outputItem);
            options.aiPanel?.addEngineMessage(outputItem);
          }
          continue;
        }

        let result = '';
        try {
          result = await executeFunction(context, functionName, language, args, options.aiPanel);
        } catch (e) {
          const msg = e instanceof Error ? e.message : String(e);
          result = `Error executing ${functionName}: ${msg}`;
        }

        const outputItem = langTool.createToolOutputMessage(callId, functionName, result);
        input.push(outputItem);
        options.aiPanel?.addEngineMessage(outputItem);

        // Clarification is a terminal action: ask the user and stop.
        if (functionName === 'ask_user_for_clarification')
          return '';
      }

      // If the model called tools, we must continue so it can consume tool outputs.
      if (hadToolCalls) {
        if (iterations >= maxIterations && options.aiPanel) {
          const cont = await options.aiPanel.addConfirmMessage('Maximum iterations reached. Do you want to continue generating the SQL query?');
          if (cont)
            maxIterations += 10; // reset max iterations to continue
        }
        continue;
      }

      const content = findLast(response.content, (c) => c.type === 'text')?.text!;
      if ((content?.length ?? 0) === 0)
        throw new Error('Model returned no text and no tool calls');

      if (content.startsWith('CLARIFICATION NEEDED:')) {
        const clarificationMsg = content.replace('CLARIFICATION NEEDED:', '').trim();
        options.aiPanel?.addUiMessage(clarificationMsg, false);
        return '';
      }

      if (content.startsWith('REPLY ONLY:')) {
        const reply = content.replace('REPLY ONLY:', '').trim();
        options.aiPanel?.addUiMessage(reply, false);
        return '';
      }
      //   let sanitizedContent = content;
      //   if (content && content)

      const validationError = validateScriptHeader(content, language);
      if (validationError) {
        const fixMsg = langTool.createUserMessage(`The generated script has an issue: ${validationError}. Please fix it and provide the corrected script. ALWAYS REMEMBER THAT YOU SHOULD OUTPUT THE CODE ONLY, WITHOUT ANY EXPLANATION, MARKDOWN, OR FORMATTING WITH CORRECT COMMENTS.`);
        input.push(fixMsg);
        options.aiPanel?.addEngineMessage(fixMsg);
        continue;
      }

      return content;
    }

    // If we exhausted iterations
    throw new Error('Failed to generate script: Maximum iterations reached or no valid script returned');
  } catch (e) {
    if (e instanceof ClarificationNeededError)
      return '';
    else if (e instanceof NoScriptGeneratedError)
      throw new Error('No script was generated based on the user request.');
    else {
      const msg = e instanceof Error ? e.message : String(e);
      throw new Error(`Script generation failed: ${msg}`);
    }
  } finally {
    abortSub.unsubscribe();
    context.dispose();
  }
}

/**
 * Generates the system prompt based on language
 */
function getSystemPrompt(language: DG.ScriptingLanguage, vectorStoreId: string | null): string {
  const languageSpecifics = getLanguageSpecifics(language);

  return `You are an expert Datagrok script developer specializing in ${language.toUpperCase()}. Your task is to generate production-ready Datagrok scripts based on user requirements.

IN CASE OF JAVASCRIPT SCRIPTS - NEVER USE TYPESCRIPT SYNTAX OR FEATURES, plain JAVASCRIPT ONLY.
YOUR RESPONSE CODE MUST BE STANDALONE, NO IMPORTS OR EXPORTS. ALL THAT IS AVAILABLE YOU CAN ASSUME IS ALREADY IMPORTED IN THE ENVIRONMENT.

Your name is Datagrok Script Expert.

CRITICAL SCRIPT FORMAT REQUIREMENTS:
${languageSpecifics.headerFormat}

AVAILABLE LIBRARIES AND APIS:
${languageSpecifics.availableLibraries}

WORKFLOW:
1. Analyze the user's request thoroughly
2. Identify required parameters and decide which must be script inputs.
  - If the user did not specify required parameters (e.g., "aggregate the table by given column"), you MUST infer that these are missing and MUST become explicit script inputs.
  - Typical examples of required parameters that should become inputs:
    * group-by column(s)
    * value/measure column(s)
    * aggregation function(s) (sum/avg/min/max/count/etc.)
    * filtering conditions, date ranges, thresholds
    * output shape (dataframe/column/graphics) when not obvious
3. If any required parameter is still ambiguous (too many plausible interpretations), ask clarifying questions and STOP.
  - Use the tool ask_user_for_clarification.
  - Ask only the minimum questions needed to proceed.
  - The clarification questions MUST directly map to future #input lines (columns/strings/numbers/booleans).
  - After asking for clarification, DO NOT generate code in the same turn.
4. Use search_documentation to understand Datagrok APIs and features (MANDATORY before using any Datagrok-specific API)
5. Use find_similar_script_samples to find code examples for similar functionality (HIGHLY RECOMMENDED)
6. (JavaScript only) Before using any DG/grok/ui member that is NOT directly present in the examples you found, verify it exists:
  - Call list_js_members on the parent object/namespace (e.g., list_js_members("grok.chem") before using grok.chem.sketcher)
  - If the member is not present, do NOT use it; pick an alternative from docs/examples or ask for clarification.
6.1 (JavaScript only, if available) Before using any Datagrok JS API method/property (e.g., grok.data.demo, DG.DataFrame.fromCsv, ui.dialog), you ALWAYS MUST (!!!!!) check it exists with search_js_api_sources.
  - This tool returns source code chunks from the datagrok js-api; use it to verify the API surface.
  - If the member does not appear in the returned source chunks, do NOT use it; pick an alternative or ask for clarification.
7. Generate the complete script with proper header comments
8. Validate the script has all required annotations
9. Return ONLY the complete script code (plain text; no markdown, no code fences, no explanations)


CRITICAL RULES:
- ALWAYS include proper header comments with #name, #description, #language, #input, and #output annotations
- Treat any user-specified or clarified parameters as #input annotations (do not hardcode them unless the user explicitly wants a fixed constant)
- NEVER hallucinate API methods or properties - ALWAYS search documentation first
- For any language: ALWAYS use find_similar_script_samples before writing code that uses DG, grok, or ui namespaces
- For JavaScript: if you use a DG/grok/ui method/property that you did not see in examples, you MUST first verify it exists using list_js_members
- Input columns should have semantic types when relevant (e.g., {semType: Molecule} for molecular data)
- Common script patterns:
  * Data transformation: input dataframe + column(s) → output dataframe
  * Column calculation: input dataframe + column(s) → output column
  * Visualization: input dataframe + params → output graphics/viewer
  * Analysis: input dataframe → output string/value
- Use precise semantic types: Molecule, Sequence, Cell Line, etc.
- ALWAYS test your understanding by searching documentation before implementation
- If you cannot generate a valid script after research, use ask_user_for_clarification${vectorStoreId ? '\n- Documentation search is MANDATORY for any Datagrok-specific functionality' : ''}${vectorStoreId && language === 'javascript' ? '\n- For JavaScript, you MUST use search_js_api_sources to verify JS API members before use' : ''}
- If the user is just asking a question or chatting (not requesting a script), respond with "REPLY ONLY: [your response]"
- Final output MUST be only the script code. Do not include markdown, bullets, backticks, headings, or any extra text.
- THE SCRIPT YOU GENERATE MUST BE a script that is run directly, not a function.
- MAKE SURE YOU DO NOT REDECLARE THE VARIABLES THAT ARE ALREADY DECLARED IN THE HEADER COMMENTS AS INPUTS (the output variable you MUST DECLARE). you can assume that input variables are already declared!!!!.
- Avalialable Variable types: dataframe, column, string, int, double, boolean, graphics, viewer, list, dynamic. DO NOT USE ANY OTHER VARIABLE TYPES (except if seen in examples)!!!!

RESPONSE MODES:
1. SCRIPT MODE (default): Return complete script code with header comments
2. REPLY ONLY MODE: Start response with "REPLY ONLY:" for non-script questions
3. CLARIFICATION MODE: Use ask_user_for_clarification tool when essential info is missing, then stop.

When you have the final script ready, respond with ONLY the complete script code including header comments (no markdown code blocks, no extra explanation).

NEVER IMPORT ANYTHING in JAVASCRIPT SCRIPTS, In python/R/Julia/Matlab scripts, use only standard libraries mentioned in available libraries section.
NEVER EXPORT ANYTHING
NEVER RETURN ANYTHING. IF YOU HAVE AN OUTPUT, You just assign your output variable as per the header comments
for example:
//output: dataframe result [Output description]
result = ...

ALWAYS USE DOCUMENTATION SEARCH AND SIMILAR SCRIPTS SEARCH BEFORE outputing final result. DO NOT ASSUME STUFF!!!!
`;
}

/**
 * Gets language-specific information
 */
function getLanguageSpecifics(language: DG.ScriptingLanguage): {
  headerFormat: string,
  availableLibraries: string,
  commentPrefix: string
} {
  switch (language) {
  case 'javascript':
    return {
      commentPrefix: '//',
      headerFormat: `JavaScript scripts MUST start with these comment annotations:
//name: [Category] | [Script Name]
//description: [Brief description]
//language: javascript
//input: dataframe table [Description of dataframe input]
//input: column col {semType: [SemanticType]} [Column description]
//output: dataframe result [Output description]

Example:
//name: Transform | Calculate Descriptors
//description: Calculates molecular descriptors
//language: javascript
//input: dataframe table
//input: column molecules {semType: Molecule}
//output: dataframe descriptors`,
      availableLibraries: `- Full Datagrok JavaScript API via DG, grok, and ui namespaces
- DG.DataFrame, DG.Column, DG.Viewer for data manipulation
- grok.shell, grok.data, grok.chem for platform features
- ui.div, ui.button, ui.input for UI components
- MANDATORY: Search documentation and find similar samples before using any API`
    };
  case 'python':
    return {
      commentPrefix: '#',
      headerFormat: `Python scripts MUST start with these comment annotations:
#name: [Category] | [Script Name]
#description: [Brief description]
#language: python
#input: dataframe table [Description]
#input: string column {semType: [Type]} [Description]
#output: dataframe result [Description]

Example:
#name: Chemistry | Gasteiger Partial Charges
#description: RDKit-based script for calculating charges
#language: python
#input: string mol = "CCO" {semType: Molecule} [Molecule in SMILES]
#output: graphics charges [The Gasteiger partial charges]`,
      availableLibraries: `- rdkit: Full RDKit library for chemistry
- pandas: DataFrame operations
- numpy: Numerical operations
- Standard Python libraries`
    };
  case 'r':
    return {
      commentPrefix: '#',
      headerFormat: `R scripts MUST start with these comment annotations:
#name: [Category] | [Script Name]
#description: [Brief description]
#language: r
#input: dataframe table
#input: column col {semType: [Type]}
#output: dataframe result`,
      availableLibraries: `- Standard R libraries and packages
- Data manipulation with data.frame, dplyr
- Statistical functions`
    };
  case 'julia':
    return {
      commentPrefix: '#',
      headerFormat: `Julia scripts MUST start with these comment annotations:
#name: [Category] | [Script Name]
#description: [Brief description]
#language: julia
#input: dataframe table
#output: dataframe result`,
      availableLibraries: `- Standard Julia packages
- DataFrames.jl for data manipulation`
    };
  case 'octave':
    return {
      commentPrefix: '%',
      headerFormat: `MATLAB scripts MUST start with these comment annotations:
%name: [Category] | [Script Name]
%description: [Brief description]
%language: matlab
%input: dataframe table
%output: dataframe result`,
      availableLibraries: `- Standard MATLAB functions and toolboxes
- Table and array operations`
    };
  default:
    throw new Error(`Unsupported language: ${language}`);
  }
}

/**
 * Defines tools available for script generation
 */
function getTools(language: DG.ScriptingLanguage, vectorStoreId: string | null) {
  const tools: LanguageModelV3FunctionTool[] = [
    {
      type: 'function',
      name: 'ask_user_for_clarification',
      description: 'Ask the user for clarification when the requirements are vague or missing critical information. Use this when you need to know: required inputs, expected outputs, specific parameters, data types, or any other essential details.',
      inputSchema: {
        type: 'object',
        properties: {
          question: {
            type: 'string',
            description: 'Clear, specific question(s) for the user in markdown format. Be professional and helpful.'
          }
        },
        required: ['question'],
        additionalProperties: false
      },
      strict: true
    },
    {
      type: 'function',
      name: 'reply_to_user',
      description: 'Send a message to the user to explain your reasoning, provide updates, or answer questions. Use this for intermediate communication, NOT for the final script.',
      inputSchema: {
        type: 'object',
        properties: {
          message: {
            type: 'string',
            description: 'Your message to the user in markdown format'
          }
        },
        required: ['message'],
        additionalProperties: false
      },
      strict: true
    }
  ];

  // Add documentation search if vector store is available
  if (vectorStoreId) {
    tools.push({
      type: 'function',
      name: 'search_documentation',
      description: 'Search Datagrok documentation and API references. MANDATORY before using any Datagrok API. Returns detailed information about classes, methods, properties, and usage examples.',
      inputSchema: {
        type: 'object',
        properties: {
          query: {
            type: 'string',
            description: 'Detailed search query describing what API or feature you need to understand. Include context about what you are trying to do.'
          },
          maxResults: {
            type: 'integer',
            description: 'Maximum number of results to return (1-50)',
            default: 4
          }
        },
        required: ['query', 'maxResults'],
        additionalProperties: false
      },
      strict: true
    });
  }

  // Add API sample search for any language
  tools.push({
    type: 'function',
    name: 'find_similar_script_samples',
    description: 'Find similar scripts code samples from Datagrok codebase. HIGHLY RECOMMENDED before writing any code using DG, grok, or ui APIs. Returns real working examples.',
    inputSchema: {
      type: 'object',
      properties: {
        description: {
          type: 'string',
          description: 'Description of the functionality you want to implement'
        }
      },
      required: ['description'],
      additionalProperties: false
    },
    strict: true
  });

  if (language === 'javascript') {
    if (vectorStoreId) {
      tools.push({
        type: 'function',
        name: 'search_js_api_sources',
        description: 'Search Datagrok JavaScript API source code (js-api) to verify that a method/property exists. Use this before calling any grok/ui/DG members that are not already confirmed by examples. Returns relevant source chunks.',
        inputSchema: {
          type: 'object',
          properties: {
            query: {
              type: 'string',
              description: 'Search query for the JS API source code (e.g., "grok.data.demo", "DG.DataFrame.fromCsv", "ui.dialog")'
            },
            maxResults: {
              type: 'integer',
              description: 'Maximum number of results to return (1-50)',
              default: 4
            }
          },
          required: ['query', 'maxResults'],
          additionalProperties: false
        },
        strict: true
      });
    }
    tools.push({
      type: 'function',
      name: 'list_js_members',
      description: 'List enumerable members (Object.keys) of a JavaScript object/namespace to verify methods/properties exist. Use this to avoid hallucinating APIs. Input must be a dotted identifier expression like "grok.chem" or "DG" or "ui".',
      inputSchema: {
        type: 'object',
        properties: {
          expression: {
            type: 'string',
            description: 'A dotted identifier expression referring to an in-scope object/namespace (e.g., "grok.chem", "grok.shell", "DG", "ui", "grok.dapi.entities", etc)'
          }
        },
        required: ['expression'],
        additionalProperties: false
      },
      strict: true
    });
  }


  return tools;
}

/**
 * Executes a function call from the model
 */
async function executeFunction(
  context: ScriptGenerationContext,
  functionName: string,
  language: DG.ScriptingLanguage,
  args: JsonObject | null,
  aiPanel?: AIPanelFuncs<MessageType>
): Promise<string> {
  switch (functionName) {
  case 'search_documentation':
    aiPanel?.addUiMessage(`Searching documentation for "${getStringProp(args, 'query') ?? ''}"...`, false);
    return await context.searchDocumentation(getStringProp(args, 'query') ?? '', getNumberProp(args, 'maxResults'));

  case 'find_similar_script_samples':
    aiPanel?.addUiMessage('Looking for similar JavaScript samples...', false);
    return await context.findSimilarScriptSamples(getStringProp(args, 'description') ?? '', language);

  case 'search_js_api_sources':
    aiPanel?.addUiMessage(`Searching JS API resources for "${getStringProp(args, 'query') ?? ''}"...`, false);
    return await context.searchJsApiSources(getStringProp(args, 'query') ?? '', getNumberProp(args, 'maxResults'));

  case 'list_js_members':
    aiPanel?.addUiMessage(`Checking JS members for "${getStringProp(args, 'expression') ?? ''}"...`, false);
    return context.listJsMembers(getStringProp(args, 'expression') ?? '');

  case 'ask_user_for_clarification':
    // Show question to user through panel
    if (aiPanel)
      aiPanel.addUiMessage(getStringProp(args, 'question') ?? '', false);

    return 'Question shown to user. Wait for user response before proceeding.';

  case 'reply_to_user':
    // Show message to user through panel
    if (aiPanel)
      aiPanel.addUiMessage(getStringProp(args, 'message') ?? '', false);

    return 'Message sent to user.';

  default:
    throw new Error(`Unknown function: ${functionName}`);
  }
}

/**
 * Validates script header format
 */
function validateScriptHeader(script: string, language: DG.ScriptingLanguage): string | null {
  const lines = script.split('\n');
  const commentPrefix = getLanguageSpecifics(language).commentPrefix;

  // Check for required annotations
  let hasName = false;
  let hasDescription = false;
  let hasLanguage = false;
  let hasOutput = false;
  let iter = 0;
  if (lines[0]?.startsWith('```'))
    return 'Script should not be enclosed in markdown code fences, there SHOULD BE NO ``` OR ANY MARKDOWN!!! CODE ONLY!!!';
  for (const line of lines.slice(0, 20)) { // Check first 20 lines
    const trimmed = line.trim();
    if (trimmed.startsWith(`${commentPrefix}name:`) && iter === 0) hasName = true;
    if (trimmed.startsWith(`${commentPrefix}description:`)) hasDescription = true;
    if (trimmed.startsWith(`${commentPrefix}language:`)) hasLanguage = true;
    if (trimmed.startsWith(`${commentPrefix}output:`)) hasOutput = true;
    iter++;
  }

  if (!hasName) return 'Missing required annotation: name. SCRIPT MUST start with name annotation.';
  if (!hasDescription) return 'Missing required annotation: description';
  if (!hasLanguage) return 'Missing required annotation: language';
  if (!hasOutput) return 'Missing required annotation: output';

  // check if last line includes return statement or output assignment
  if (lines[lines.length - 1].includes('return '))
    return 'Script should not use return statement. Assign output variable as per header comments.';
  // Input is optional for some scripts

  return null; // Valid
}

/**
 * Context class that manages tool execution for script generation
 */
class ScriptGenerationContext {
  constructor(private vectorStoreId: string | null) {}

  dispose() {
    // Cleanup if needed
  }

  async searchDocumentation(query: string, maxResults: number | null, filterOptions?: {
    fileExtension?: string,
    firstFolder?: string, // useful for specific searches, like "js-api", "libraries", "packages" etc.
    secondFolder?: string, // useful for specific searches, inside of the firstFolder, like specific plugin, library, etc.
  }): Promise<string> {
    if (!this.vectorStoreId)
      return 'Documentation search is not available. Vector store ID not configured.';

    if (!query || query.trim().length === 0)
      return 'Error: empty documentation search query.';

    try {
      const openai = LLMClient.getInstance().openai;
      const filterObject = {filters: [] as ComparisonFilter[], type: 'and'} satisfies CompoundFilter;
      if (filterOptions) {
        if (filterOptions.fileExtension) {
          filterObject.filters.push({
            key: 'originalExtension',
            value: filterOptions.fileExtension,
            type: 'eq',
          });
        }
        if (filterOptions.firstFolder) {
          filterObject.filters.push({
            key: 'firstParentFolder',
            value: filterOptions.firstFolder,
            type: 'eq',
          });
        }
        if (filterOptions.secondFolder) {
          filterObject.filters.push({
            key: 'secondParentFolder',
            value: filterOptions.secondFolder,
            type: 'eq',
          });
        }
      }
      const searchData = await openai.vectorStores.search(this.vectorStoreId, {
        query,
        max_num_results: maxResults ?? 4,
        rewrite_query: false,
        filters: filterObject.filters.length > 0 ? filterObject : undefined,
      });

      if (!searchData.data || searchData.data.length === 0)
        return 'No relevant documentation found for this query.';

      let formattedResults = `Found ${searchData.data.length} relevant documentation sections:\n\n`;
      for (const item of searchData.data) {
        formattedResults += `**Source**: ${item.filename} (score: ${item.score.toFixed(2)})\n`;
        for (const content of item.content) {
          if (content.type === 'text')
            formattedResults += `${content.text}\n\n`;
        }
        formattedResults += '---\n\n';
      }
      return formattedResults;
    } catch (e) {
      const msg = e instanceof Error ? e.message : String(e);
      return `Error searching documentation: ${msg}`;
    }
  }

  async findSimilarScriptSamples(description: string, language: DG.ScriptingLanguage): Promise<string> {
    const similarSampels = await grok.ai.searchEntities(description, 0.3, 10, ['Script']);
    if (similarSampels.length === 0)
      return 'No similar script samples found.';
    const scripts = similarSampels
      .filter((e) => e instanceof DG.Script && e.language === language) as DG.Script[];
    if (scripts.length === 0)
      return 'No similar script samples found.';
    let result = `Found ${scripts.length} similar script samples:\n\n`;
    for (const script of scripts) {
      result += `**Name**: ${script.name}\n`;
      result += `**Description**: ${script.description}\n`;
      result += `**Code**:\n${script.script}\n\n---\n\n`;
    }
    return result;
  }

  async searchJsApiSources(query: string, maxResults: number | null): Promise<string> {
    if (!this.vectorStoreId)
      return 'JS API source search is not available. Vector store ID not configured.';

    if (!query || query.trim().length === 0)
      return 'Error: empty JS API source search query.';

    return await this.searchDocumentation(query, maxResults, {
      firstFolder: 'js-api',
    });
  }

  listJsMembers(expression: string): string {
    if (!expression || expression.trim().length === 0)
      return 'Error: expression is required.';
    const expr = expression.trim();
    if (!isSafeJsMemberExpression(expr))
      return 'Error: unsafe expression. Only dotted identifiers like "grok.chem" are allowed. consider use js-api search instead';

    try {
      // eslint-disable-next-line no-eval
      const value = eval(expr) as object | null | undefined;
      if (value === null || value === undefined)
        return 'Error: expression evaluated to null/undefined.';

      const keys = Object.keys(value as Record<string, string | number | boolean | null | object>);
      return JSON.stringify(keys);
    } catch (e) {
      const msg = e instanceof Error ? e.message : String(e);
      return `Error evaluating expression: ${msg}`;
    }
  }
}
