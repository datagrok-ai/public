/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelOption, ModelType, LLMClient} from './LLM-client';
import {ChatModel} from 'openai/resources/shared';
import {findLast, getAIAbortSubscription} from '../utils';
import {AIPanelFuncs, MessageType} from './panel';
import {LLMCredsManager} from './creds';
import {LanguageModelV3, LanguageModelV3FunctionTool, LanguageModelV3Message} from '@ai-sdk/provider';

type JsonPrimitive = string | number | boolean | null;
type JsonValue = JsonPrimitive | JsonObject | JsonValue[];
type JsonObject = { [key: string]: JsonValue };

type ViewerProps = Record<string, JsonValue>;

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

/**
 * Handles AI-assisted table view interactions with ask and agent modes
 * @throws Error if processing fails
 * @param prompt - User's natural language query
 * @param tableView - The table view to work with
 * @param mode - 'ask' for highlighting elements or 'agent' for full viewer manipulation
 * @param options - Additional options including message history and AI panel functions
 */
export async function processTableViewAIRequest(
  prompt: string,
  tableView: DG.TableView,
  options: {
    oldMessages?: MessageType[],
    aiPanel?: AIPanelFuncs<MessageType>,
    modelType?: ModelOption,
  } = {}
): Promise<void> {
  let aborted = false;
  // Create tool execution context
  const context = new TableViewContext(tableView);
  const abortSub = getAIAbortSubscription().subscribe(() => {
    aborted = true;
    console.log('Aborting table view AI processing as per user request');
    abortSub.unsubscribe();
  });

  try {
    // Initialize OpenAI client
    const langTool = LLMClient.getInstance();
    const modelType = options.modelType ?? 'Fast';
    const client = langTool.aiModels[modelType];
    // Always use agent mode workflow with conditional system prompt
    await processAgentMode(client, langTool, prompt, context, options, () => aborted);
  } finally {
    abortSub.unsubscribe();
  }
}

/**
 * Unified agent workflow with conditional system prompt based on mode
 */
async function processAgentMode(
  client: LanguageModelV3,
  langTool: LLMClient,
  prompt: string,
  context: TableViewContext,
  options: {
    oldMessages?: MessageType[],
    aiPanel?: AIPanelFuncs<MessageType>,
    modelType?: ModelOption
  },
  isAborted: () => boolean,
): Promise<void> {
  let maxIterations = 15;
  let iterations = 0;

  const vectorStoreId = LLMCredsManager.getVectorStoreId();

  // Conditional system prompt based on mode
  const systemMessage =
    `You are Datagrok Table View Assistant, an expert in data visualization and analysis using Datagrok viewers.

Your name is Datagrok Table View Assistant.

WORKFLOW:
1. You have access to the current table structure and available viewer types
2. Use describe_viewer to understand viewer properties before adding/modifying
3. Use add_viewer to create new visualizations
4. Use adjust_viewer to modify existing viewers (if you know the ID) or find_viewers_by_type to find viewer IDs
5. Use list_current_viewers to see what's currently displayed
6. Iterate and refine based on user feedback

CRITICAL RULES:
- ALWAYS use describe_viewer before adding a viewer to understand its properties
- When setting viewer properties, use exact column names from the table structure
- For categorical columns, leverage the category information provided
- Pay attention to semantic types (semType) and units when relevant
- Use list_current_viewers to understand the current state before making changes
- If user mentions "this viewer" or "that chart", use find_viewers_by_type to identify which one
- If adjust_viewer fails because viewer doesn't exist, add it instead
- Format your responses in markdown for better readability
- If the user asks a general question that doesn't require viewer manipulation, just answer it conversationally
- When user request is ambiguous, ask for clarification
- ALWAYS use highlight_element to guide users to UI elements when relevant!${vectorStoreId ? '\n- DOCUMENTATION SEARCH: Use search_documentation tool IMMEDIATELY when the user asks about any Datagrok features, functionality, or workflows that are not directly related to adding/modifying viewers in the current table. This includes questions like "how do I do X" or "where is feature Y". NEVER invent or assume information about Datagrok features, menu locations, package names, or workflows - ALWAYS search documentation first! If documentation search returns no results, tell the user you cannot find information rather than making assumptions.' : ''}
- STRICT NO-INVENTION POLICY: Do NOT make up menu locations, feature names, dialog boxes, keyboard shortcuts, package names, or step-by-step instructions unless they come directly from documentation search results or are about the specific table/viewer tools you have access to.

Your responses should be informative, explaining what you're doing and why.`;

  const input: MessageType[] = options.oldMessages ? [...options.oldMessages] : [];

  if (input.length === 0) {
    const tableDescription = context.getTableDescription();
    const availableViewers = context.getAvailableViewers();

    const initialContext = `${tableDescription}\n\n${availableViewers}`;

    const systemMsg = langTool.createSystemMessage(systemMessage);
    const userMsg = langTool.createUserMessage(`${initialContext}\n\nUser request: ${prompt}`);
    input.push(systemMsg);
    input.push(userMsg);
    options.aiPanel?.addEngineMessage(systemMsg);
    options.aiPanel?.addUserMessage(userMsg, prompt);
  } else {
    const userMsg = langTool.createUserMessage(prompt);
    input.push(userMsg);
    options.aiPanel?.addUserMessage(userMsg, prompt);
  }

  const tools: LanguageModelV3FunctionTool[] = [
    {
      type: 'function',
      name: 'list_ui_elements',
      description: 'Get a list of all available UI elements in the table view with their IDs and descriptions. Use this when users ask where to find features or how to access functionality. Call this before highlight_element to know which elements are available.',
      inputSchema: {
        type: 'object',
        properties: {},
        required: [],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'highlight_element',
      description: 'Highlights a UI element to guide the user. Use this when you want to show users where specific features or controls are located. Use list_ui_elements first to get valid element ID or names.',
      inputSchema: {
        type: 'object',
        properties: {
          elementId: {
            type: 'string',
            description: 'ID or name of the UI element to highlight'
          },
          description: {
            type: 'string',
            description: 'Brief description of why this element is being highlighted'
          }
        },
        required: ['elementId', 'description'],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'get_available_viewers',
      description: 'Get a list of all available viewer types that can be added to the table view. Use this if you forget which viewer types are available or need to remind yourself of the options.',
      inputSchema: {
        type: 'object',
        properties: {},
        required: [],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'describe_viewer',
      description: 'Get detailed information about a viewer type including all its properties and their types. ALWAYS call this BEFORE adding a new viewer to understand what properties are available and their correct names/types. This describes the viewer type schema, not existing viewer instances.',
      inputSchema: {
        type: 'object',
        properties: {
          viewerType: {
            type: 'string',
            description: 'Type of the viewer to describe (e.g., "Histogram", "Scatter plot"). must match exactly one of the available viewer types supplied.'
          }
        },
        required: ['viewerType'],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'add_viewer',
      description: 'Add a new viewer to the table view with specified properties. IMPORTANT: Set all desired properties during creation using the viewerProperties parameter - avoid creating a viewer with empty properties and then adjusting it separately. Use describe_viewer first to understand available properties.',
      inputSchema: {
        type: 'object',
        properties: {
          viewerType: {
            type: 'string',
            description: 'Type of viewer to add'
          },
          viewerProperties: {
            type: 'object',
            description: 'Properties to configure for the viewer (e.g., {"valueColumnName": "Age", "markerSize": 15})',
            properties: {},
            required: [],
            patternProperties: {
              '.*': {type: 'string'},
            },
            additionalProperties: false,
          }
        },
        required: ['viewerType', 'viewerProperties'],
        additionalProperties: false,
      },
      // strict: true,
    },
    {
      type: 'function',
      name: 'find_viewers_by_type',
      description: 'Find all viewers of a specific type that are already added to the table view and get their IDs and current properties. Use this when you need to identify existing viewers of a certain type before modifying them, or when the user refers to "the histogram" or "the scatter plot".',
      inputSchema: {
        type: 'object',
        properties: {
          viewerType: {
            type: 'string',
            description: 'Type of viewers to find'
          }
        },
        required: ['viewerType'],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'adjust_viewer',
      description: 'Adjust properties of an EXISTING viewer by its ID. Use this to modify viewers that are already created, not for initial configuration. If viewer is not found, will inform you to add it instead. For new viewers, use add_viewer with all properties set from the start.',
      inputSchema: {
        type: 'object',
        properties: {
          viewerId: {
            type: 'string',
            description: 'ID of the viewer to adjust'
          },
          viewerProperties: {
            type: 'object',
            description: 'Properties to update',
            properties: {},
            required: [],
            patternProperties: {
              '.*': {type: 'string'},
            },
            additionalProperties: false,
          }
        },
        required: ['viewerId', 'viewerProperties'],
        additionalProperties: false,
      },
      // strict: true,
    },
    {
      type: 'function',
      name: 'list_current_viewers',
      description: 'Get a list of all currently added viewers with their IDs, types, and configured properties. Use this to see what viewers already exist before adding new ones or to find viewer IDs for adjustment.',
      inputSchema: {
        type: 'object',
        properties: {},
        required: [],
        additionalProperties: false,
      },
      strict: true,
    },
    {
      type: 'function',
      name: 'reply_to_user',
      description: 'Send a message to the user to explain reasoning, ask for clarification, or provide updates',
      inputSchema: {
        type: 'object',
        properties: {
          message: {
            type: 'string',
            description: 'Message to send to the user'
          }
        },
        required: ['message'],
        additionalProperties: false,
      },
      strict: true,
    }
  ];

  // Add search_documentation tool only if vectorStoreId is provided
  if (vectorStoreId) {
    tools.push({
      type: 'function',
      name: 'search_documentation',
      description: 'Search external Datagrok documentation for authoritative information. CRITICAL: Use this tool IMMEDIATELY when user asks ANY question about: 1) How to access Datagrok features/functionality, 2) Where to find menu items or commands, 3) How to perform operations not directly related to viewer manipulation, 4) Package names or requirements, 5) Keyboard shortcuts, 6) Dialogs or workflows. DO NOT make assumptions or invent UI locations, feature names, or steps - ALWAYS search documentation first! If you are even slightly uncertain about something, USE THIS TOOL. Only skip this tool for questions that are purely about manipulating viewers in the current table.',
      inputSchema: {
        type: 'object',
        properties: {
          query: {
            type: 'string',
            description: 'The search query for finding relevant documentation'
          },
          maxResults: {
            type: 'integer',
            description: 'Maximum number of results to return (1-50)',
            default: 5
          }
        },
        required: ['query', 'maxResults'],
        additionalProperties: false,
      },
      strict: true,
    });
  }

  while (iterations < maxIterations) {
    if (isAborted()) {
      options.aiPanel?.addUiMessage('**Processing aborted by user**', false);
      return;
    }

    iterations++;
    const response = await client.doGenerate({
      prompt: input,
      tools,
      providerOptions: {
        openai: {
          ...(ModelType.Coding.startsWith('gpt-5') ? {reasoning: {effort: 'medium'}} : {}),
        }
      }
    });

    // The AI SDK's doGenerate puts item IDs in providerMetadata, but the
    // Responses API input converter only reads reasoning items from providerOptions
    // (tool-call has a fallback to providerMetadata, but reasoning does not).
    // Copy providerMetadata → providerOptions so reasoning item_references are created.
    const fixedContent = response.content.map((item) => {
      const patched = {
        ...item,
        providerOptions: (item as any).providerOptions ?? (item as any).providerMetadata,
      };
      if (item.type === 'tool-call') {
        return {
          ...patched,
          input: typeof item.input === 'string' ? parseJsonObject(item.input) ?? item.input : item.input
        };
      }
      return patched;
    });

    const formattedOutput: LanguageModelV3Message = {
      role: 'assistant',
      // @ts-ignore
      content: fixedContent
    };
    //const outputs: MessageType[] = response.content.map((item) => ({role: 'assistant', content: [item]} as MessageType));
    input.push(formattedOutput);
    options.aiPanel?.addEngineMessage(formattedOutput);

    let hadToolCalls = false;
    for (const call of response.content.filter((c) => c.type === 'tool-call')) {
      hadToolCalls = true;

      const functionName = call.toolName;
      const callId = call.toolCallId;
      const rawArgs = call.input ?? '{}';
      const argsObj = parseJsonObject(rawArgs);

      if (!functionName || !callId) {
        const outputItem = langTool.createToolOutputMessage(callId ?? '', functionName ?? '', 'Error: malformed function call (missing name/call_id).');
        if (callId) {
          input.push(outputItem);
          options.aiPanel?.addEngineMessage(outputItem);
        }
        continue;
      }

      let result: string;
      try {
        switch (functionName) {
        case 'list_ui_elements':
          result = await context.listUIElements();
          options.aiPanel?.addUiMessage('🔍 Looking for available UI elements', false);
          break;
        case 'get_available_viewers':
          result = context.getAvailableViewers();
          options.aiPanel?.addUiMessage('📊 Getting available viewer types', false);
          break;
        case 'highlight_element': {
          result = 'Element highlighted';
          const elementId = getStringProp(argsObj, 'elementId') ?? '';
          const description = getStringProp(argsObj, 'description') ?? '';
          const actRes = await context.highlightElement(elementId, description);
          options.aiPanel?.addUiMessage(`Highlighting: **${elementId}**`, false);
          if (actRes != null) {
            // if it is not null, then it is an error message
            result = actRes;
            options.aiPanel?.addUiMessage(`⚠️ ${actRes}`, false);
          }
          break;
        }
        case 'describe_viewer': {
          const viewerType = getStringProp(argsObj, 'viewerType') ?? '';
          result = await context.describeViewer(viewerType);
          options.aiPanel?.addUiMessage(`📋 Getting available properties for **${viewerType}** viewer`, false);
          break;
        }
        case 'add_viewer': {
          const viewerType = getStringProp(argsObj, 'viewerType') ?? '';
          const viewerProps = (isJsonObject(argsObj?.viewerProperties) ? argsObj.viewerProperties :
            isJsonObject(argsObj?.properties) ? argsObj.properties : {}) as ViewerProps;
          result = await context.addViewer(viewerType, viewerProps);
          options.aiPanel?.addUiMessage(`➕ Added **${viewerType}** viewer ${viewerProps && Object.keys(viewerProps).length > 0 ? `with properties: \n\`\`\`json\n${JSON.stringify(viewerProps, null, 2)}\n\`\`\`` : ''}`, false);
          break;
        }
        case 'find_viewers_by_type': {
          const viewerType = getStringProp(argsObj, 'viewerType') ?? '';
          options.aiPanel?.addUiMessage(`🔎 Searching for **${viewerType}** viewers`, false);
          result = await context.findViewersByType(viewerType);
          break;
        }
        case 'adjust_viewer': {
          const viewerId = getStringProp(argsObj, 'viewerId') ?? '';
          const viewerProps = (isJsonObject(argsObj?.viewerProperties) ? argsObj.viewerProperties :
            isJsonObject(argsObj?.properties) ? argsObj.properties : {}) as ViewerProps;
          result = await context.adjustViewer(viewerId, viewerProps);
          options.aiPanel?.addUiMessage(`⚙️ Setting properties to the viewer: \n\`\`\`json\n${JSON.stringify(viewerProps, null, 2)}\n\`\`\``, false);
          break;
        }
        case 'list_current_viewers':
          result = await context.listCurrentViewers();
          options.aiPanel?.addUiMessage('📝 Listing currently open viewers', false);
          break;
        case 'reply_to_user': {
          const message = getStringProp(argsObj, 'message') ?? '';
          result = 'Message sent to user';
          options.aiPanel?.addUiMessage(`💬 ${message}`, false);
          break;
        }
        case 'search_documentation': {
          const query = getStringProp(argsObj, 'query') ?? '';
          const maxResults = getNumberProp(argsObj, 'maxResults') ?? 4;
          // Ask for user confirmation before searching external docs (example of using confirmations)
          // const confirmSearch = options.aiPanel ? await options.aiPanel.addConfirmMessage('Datagrok wants to search platform documentation for additional context. Do you want to allow this action?') : true;
          // if (!confirmSearch) {
          //   result = 'Documentation search prohibited by user. Please answer based on the available context without using external documentation.';
          //   options.aiPanel?.addUiMessage('⚠️ External documentation search prohibited by user.', false);
          //   break;
          // }
          options.aiPanel?.addUiMessage(`🔍 Searching documentation for: **${query}**`, false);
          try {
            const searchData = await langTool.openai.vectorStores.search(vectorStoreId, {
              query,
              max_num_results: maxResults,
              rewrite_query: false,
            });

            // Format the search results
            if (searchData.data && searchData.data.length > 0) {
              let formattedResults = `Found ${searchData.data.length} relevant documentation sections:\n\n`;
              for (const item of searchData.data) {
                formattedResults += `**Source**: ${item.filename} (score: ${item.score.toFixed(2)})\n`;
                for (const content of item.content) {
                  if (content.type === 'text')
                    formattedResults += `${content.text}\n\n`;
                }
                formattedResults += '---\n\n';
              }
              result = formattedResults;
              options.aiPanel?.addUiMessage(`✅ Found ${searchData.data.length} relevant documentation sections`, false);
            } else {
              result = 'No relevant documentation found for this query.';
              options.aiPanel?.addUiMessage('ℹ️ No relevant documentation found', false);
            }
          } catch (error) {
            const msg = error instanceof Error ? error.message : String(error);
            result = `Error searching documentation: ${msg}`;
            options.aiPanel?.addUiMessage(`⚠️ Error searching documentation: ${msg}`, false);
          }
          break;
        }
        default:
          result = `Unknown function: ${functionName}`;
        }
      } catch (error) {
        const msg = error instanceof Error ? error.message : String(error);
        result = `Error: ${msg}`;
        options.aiPanel?.addUiMessage(`⚠️ Error executing ${functionName}: ${msg}`, false);
      }

      const outputItem = langTool.createToolOutputMessage(callId, functionName, result);
      input.push(outputItem);
      options.aiPanel?.addEngineMessage(outputItem);
    }

    if (hadToolCalls) {
      if (iterations >= maxIterations) {
        const continueProcessing = await new Promise<boolean>((resolve) => {
          ui.dialog('Maximum Iterations Reached')
            .add(ui.divText(`The AI has reached the maximum iteration limit (${maxIterations}). Would you like to continue for 10 more iterations?`))
            .onOK(() => resolve(true))
            .onCancel(() => resolve(false))
            .show();
        });

        if (continueProcessing) {
          maxIterations += 10;
          options.aiPanel?.addUiMessage(`ℹ️ Continuing processing for ${maxIterations - iterations} more iterations...`, false);
        } else {
          options.aiPanel?.addUiMessage('⚠️ Processing stopped by user.', false);
          break;
        }
      }
      continue;
    }

    const content = findLast(response.content, (c) => c.type === 'text')?.text;
    if ((content?.length ?? 0) > 0) {
      options.aiPanel?.addUiMessage(content!, false, {
        finalResult: content
      });
      break;
    }
  }

  if (iterations >= maxIterations)
    options.aiPanel?.addUiMessage('⚠️ Maximum iteration limit reached. Please refine your request.', false);
}

const AI_GENERATED_IDENTIFIER_TAG = 'AI-generated-identifier';

/**
 * Context class that manages tool execution for table view operations
 */
class TableViewContext {
  private dummy: DG.DataFrame;
  constructor(private tableView: DG.TableView) {
    // Initialize with existing viewers
    this.dummy = grok.data.testData('demog', 100);
  }

  getTableDescription(): string {
    const df = this.tableView.dataFrame;
    let description = `## Current Table: ${df.name}\n\n`;
    description += `**${df.columns.length} Columns:**\n`;

    for (const col of df.columns.toList()) {
      description += `- **${col.name}** (${col.type})`;

      if (col.semType)
        description += ` - semType: ${col.semType}`;

      if (col.meta.units)
        description += ` - units: ${col.meta.units}`;

      if (col.isCategorical)
        description += ` - categorical (${col.categories.length} categories)`;

      description += '\n';
    }

    return description;
  }

  getAvailableViewers(): string {
    let description = `## Available Viewer Types:\n\n`;
    description += DG.Viewer.CORE_VIEWER_TYPES.join(', ');
    description += '\n\nUse describe_viewer to get detailed information about any viewer type.';
    return description;
  }

  async describeViewer(viewerType: string): Promise<string> {
    const commonProps = ['rowSource', 'filter', 'table'];
    try {
      const viewer = DG.Viewer.fromType(viewerType, this.dummy);
      let result = `Viewer properties for ${viewerType}:\n\n`;
      for (const prop of viewer.getProperties()) {
        if (!commonProps.includes(prop.name))
          result += `${prop.name} - ${prop.type} - (${prop.description || 'No Description'})` + '\n';
      }
      return result;
    } catch (e) {
      console.error(e);
      return `Error describing viewer ${viewerType}: ${e}`;
    }
  }

  getViewerIdentifier(viewer: DG.Viewer): string {
    // use viewer tags to identify viewers
    const getRndHash = () => Math.random().toString(36).substring(2, 8);
    if (!viewer.tags.get(AI_GENERATED_IDENTIFIER_TAG)) {
      // generate a random and unique identifier
      const identifier = `viewer-${getRndHash()}-${getRndHash()}-${getRndHash()}`;
      viewer.tags.set(AI_GENERATED_IDENTIFIER_TAG, identifier);
    }
    return viewer.tags.get(AI_GENERATED_IDENTIFIER_TAG)!;
  }

  async addViewer(viewerType: string, properties: ViewerProps): Promise<string> {
    const newViewer = this.tableView.addViewer(viewerType, properties as Record<string, string | number | boolean | null | object>);
    const id = this.getViewerIdentifier(newViewer);
    return `Added viewer with ID: ${id}`;
  }

  async findViewersByType(viewerType: string): Promise<string> {
    const foundViewer = Array.from(this.tableView.viewers).filter((v) => v.type === viewerType);
    if (foundViewer.length === 0)
      return `No viewers of type ${viewerType} found.`;
    const resultLines = foundViewer.map((v) => {
      const id = this.getViewerIdentifier(v);
      const props = v.getOptions().look;
      return `- ID: ${id} \n Configured Properties: \n\`\`\`json\n${JSON.stringify(props, null, 2)}\n\`\`\``;
    });
    return `Found viewers of type ${viewerType}: \n\n ${resultLines.join('\n\n')}`;
  }

  async adjustViewer(viewerId: string, properties: ViewerProps): Promise<string> {
    const foundViewer = Array.from(this.tableView.viewers).find((v) => this.getViewerIdentifier(v) === viewerId);
    if (!foundViewer)
      return `Viewer with ID ${viewerId} not found. Please add a new viewer instead.`;
    foundViewer.setOptions(properties as Record<string, string | number | boolean | null | object>);
    return `Adjusted viewer ${viewerId}`;
  }

  async listCurrentViewers(): Promise<string> {
    const allViewers = Array.from(this.tableView.viewers);
    if (allViewers.length === 0)
      return `No viewers currently added.`;
    const lines = allViewers.map((v) => {
      const id = this.getViewerIdentifier(v);
      const type = v.type;
      const props = v.getOptions().look;
      return `- ID: ${id}, Type: ${type}\n Configured Properties:\n\`\`\`json\n${JSON.stringify(props, null, 2)}\n\`\`\``;
    });
    return `Currently added viewers:\n\n${lines.join('\n\n')}`;
  }

  private highlightableElements = [
    {id: 'Toggle filters button', description: 'Button to toggle the filter panel visibility', selector: 'div.d4-ribbon .grok-icon.grok-icon-filter'},
    {id: 'Add viewer button', description: 'Button that lets you add a new viewer to the table and has a list of all viewers', selector: 'div.d4-ribbon .grok-icon.svg-icon.svg-add-viewer'},
    {id: 'Download button', description: 'Button to download the current table in various formats (csv, excel, json, parquet, sdf, etc.)', selector: 'div.d4-ribbon .grok-icon.fal.fa-arrow-to-bottom'},
    {id: 'Add new column button', description: 'Button to add a new column to the current table using various methods (from formula, from mapping, constant etc.)', selector: 'div.d4-ribbon .grok-icon.svg-icon.svg-add-new-column'},
  ];

  private getUiElements() {
    const uiElements: {id: string, description: string, element: HTMLElement}[] = this.highlightableElements
      .map((el) => ({element: document.querySelector(el.selector) as HTMLElement, id: el.id, description: el.description}))
      .filter((el) => el.element != null);
    // also add top menu items
    ['Select', 'Edit', 'View', 'Data', 'ML', 'Bio', 'Chem'].forEach((menuText) => {
      const menuItem = Array.from(this.tableView.ribbonMenu.root.querySelectorAll(`div.d4-menu-item.d4-menu-group.d4-menu-item-horz`))
        .find((el) => (el as HTMLDivElement).innerText === menuText);
      if (menuItem)
        uiElements.push({id: `${menuText} top menu`, description: `Top menu item: "${menuText}" containing relevant functions `, element: menuItem as HTMLElement});
    });
    return uiElements;
  }

  async listUIElements(): Promise<string> {
    // TODO: Populate with actual UI elements from the table view
    const elements = this.getUiElements();

    let result = 'Available UI elements:\n\n';
    for (const elem of elements)
      result += `- **${elem.id}**: ${elem.description}\n`;

    return result;
  }

  async highlightElement(elementId: string, _description: string): Promise<string | null> {
    const element = this.getUiElements().find((el) => el.id.toLowerCase() === elementId.toLowerCase());
    if (!element)
      return `Element with ID ${elementId} not found or cannot be highlighted.`;
    const domElement = element.element;
    if (!domElement || !document.body.contains(domElement))
      return `Element with ID ${elementId} not found in the current UI.`;
    ui.hints.addHintIndicator(domElement, true, 15000);
    return null;
  }
}
