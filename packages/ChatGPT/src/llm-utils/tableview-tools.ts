/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenAI} from 'openai';
import {OpenAIClient} from './openAI-client';
import {ChatModel} from 'openai/resources/shared';
import {getAIAbortSubscription} from '../utils';
import {UIMessageOptions} from './panel';

type AIPanelFuncs = {
  addUserMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam, msg: string) => void,
  addAIMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam, title: string, msg: string) => void,
  addEngineMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam) => void, // one that is not shown in the UI
  addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => void
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
  mode: 'ask' | 'agent',
  options: {
    oldMessages?: OpenAI.Chat.ChatCompletionMessageParam[]
    aiPanel?: AIPanelFuncs,
    modelName?: ChatModel
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
    const openai = OpenAIClient.getInstance().openai;
    // Always use agent mode workflow with conditional system prompt
    await processAgentMode(openai, prompt, context, options, () => aborted, mode);
  } finally {
    abortSub.unsubscribe();
  }
}

/**
 * Unified agent workflow with conditional system prompt based on mode
 */
async function processAgentMode(
  openai: OpenAI,
  prompt: string,
  context: TableViewContext,
  options: {
    oldMessages?: OpenAI.Chat.ChatCompletionMessageParam[]
    aiPanel?: AIPanelFuncs,
    modelName?: ChatModel
  },
  isAborted: () => boolean,
  mode: 'ask' | 'agent' = 'agent'
): Promise<void> {
  let maxIterations = 15;
  let iterations = 0;

  // Conditional system prompt based on mode
  const systemMessage = mode === 'ask' ?
    `You are Datagrok Table View Assistant, an expert in helping users navigate and understand the Datagrok table view interface.

Your name is Datagrok Table View Assistant.

PREFERRED APPROACH:
- Focus on providing clear, direct answers to user questions
- When users ask "where is X" or "how do I find Y", use list_ui_elements to see available UI elements, then use highlight_element to guide them
- Keep responses concise and helpful
- Prefer answering questions directly over complex viewer manipulations
- Use reply_to_user to provide explanations and guidance

AVAILABLE TOOLS:
- list_ui_elements: Get a list of all available UI elements with their IDs and descriptions (use this first when helping users find features)
- highlight_element: Show users where specific features or controls are located (use after identifying the element from list_ui_elements)
- reply_to_user: Provide clear explanations and answers
- Other viewer tools are available but should be used sparingly in ask mode

Your responses should be clear, concise, and user-friendly.` :
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
- ALWAYS use highlight_element to guide users to UI elements when relevant!

Your responses should be informative, explaining what you're doing and why.`;

  const messages: OpenAI.Chat.ChatCompletionMessageParam[] = options.oldMessages ? options.oldMessages.slice() : [];

  if (messages.length === 0) {
    const tableDescription = context.getTableDescription();
    const availableViewers = context.getAvailableViewers();

    const initialContext = `${tableDescription}\n\n${availableViewers}`;

    messages.push({role: 'system', content: systemMessage});
    messages.push({role: 'user', content: `${initialContext}\n\nUser request: ${prompt}`});

    options.aiPanel?.addUserMessage({role: 'user', content: prompt}, prompt);
  } else {
    messages.push({role: 'user', content: prompt});
    options.aiPanel?.addUserMessage({role: 'user', content: prompt}, prompt);
  }

  const tools: OpenAI.Chat.ChatCompletionTool[] = [
    {
      type: 'function',
      function: {
        name: 'list_ui_elements',
        description: 'Get a list of all available UI elements in the table view with their IDs and descriptions. Use this when users ask where to find features or how to access functionality. Call this before highlight_element to know which elements are available.',
        parameters: {
          type: 'object',
          properties: {}
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'highlight_element',
        description: 'Highlights a UI element to guide the user. Use this when you want to show users where specific features or controls are located. Use list_ui_elements first to get valid element ID or names.',
        parameters: {
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
          required: ['elementId', 'description']
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'get_available_viewers',
        description: 'Get a list of all available viewer types that can be added to the table view. Use this if you forget which viewer types are available or need to remind yourself of the options.',
        parameters: {
          type: 'object',
          properties: {}
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'describe_viewer',
        description: 'Get detailed information about a viewer type including all its properties and their types. ALWAYS call this BEFORE adding a new viewer to understand what properties are available and their correct names/types. This describes the viewer type schema, not existing viewer instances.',
        parameters: {
          type: 'object',
          properties: {
            viewerType: {
              type: 'string',
              description: 'Type of the viewer to describe (e.g., "Histogram", "Scatter plot"). must match exactly one of the available viewer types supplied.'
            }
          },
          required: ['viewerType']
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'add_viewer',
        description: 'Add a new viewer to the table view with specified properties. IMPORTANT: Set all desired properties during creation using the viewerProperties parameter - avoid creating a viewer with empty properties and then adjusting it separately. Use describe_viewer first to understand available properties.',
        parameters: {
          type: 'object',
          properties: {
            viewerType: {
              type: 'string',
              description: 'Type of viewer to add'
            },
            viewerProperties: {
              type: 'object',
              description: 'Properties to configure for the viewer (e.g., {"valueColumnName": "Age", "markerSize": 15})'
            }
          },
          required: ['viewerType', 'viewerProperties']
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'find_viewers_by_type',
        description: 'Find all viewers of a specific type that are already added to the table view and get their IDs and current properties. Use this when you need to identify existing viewers of a certain type before modifying them, or when the user refers to "the histogram" or "the scatter plot".',
        parameters: {
          type: 'object',
          properties: {
            viewerType: {
              type: 'string',
              description: 'Type of viewers to find'
            }
          },
          required: ['viewerType']
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'adjust_viewer',
        description: 'Adjust properties of an EXISTING viewer by its ID. Use this to modify viewers that are already created, not for initial configuration. If viewer is not found, will inform you to add it instead. For new viewers, use add_viewer with all properties set from the start.',
        parameters: {
          type: 'object',
          properties: {
            viewerId: {
              type: 'string',
              description: 'ID of the viewer to adjust'
            },
            viewerProperties: {
              type: 'object',
              description: 'Properties to update'
            }
          },
          required: ['viewerId', 'viewerProperties']
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'list_current_viewers',
        description: 'Get a list of all currently added viewers with their IDs, types, and configured properties. Use this to see what viewers already exist before adding new ones or to find viewer IDs for adjustment.',
        parameters: {
          type: 'object',
          properties: {}
        }
      }
    },
    {
      type: 'function',
      function: {
        name: 'reply_to_user',
        description: 'Send a message to the user to explain reasoning, ask for clarification, or provide updates',
        parameters: {
          type: 'object',
          properties: {
            message: {
              type: 'string',
              description: 'Message to send to the user'
            }
          },
          required: ['message']
        }
      }
    }
  ];

  while (iterations < maxIterations) {
    if (isAborted()) {
      options.aiPanel?.addUiMessage('**Processing aborted by user**', false);
      return;
    }

    iterations++;

    const response = await openai.chat.completions.create({
      model: options?.modelName ?? 'gpt-4o-mini',
      temperature: 1,
      reasoning_effort: options?.modelName?.startsWith('gpt-5') ? 'high' : undefined,
      messages,
      tools,
    });

    const assistantMessage = response.choices[0].message;
    messages.push(assistantMessage);
    options.aiPanel?.addEngineMessage(assistantMessage);
    // Handle text response
    if (assistantMessage.content && assistantMessage.content.trim().length > 0) {
      options.aiPanel?.addUiMessage(assistantMessage.content, false, {
        result: {finalResult: assistantMessage.content}
      });
    }

    // Handle tool calls
    const toolCalls = assistantMessage.tool_calls?.filter((tc) => tc.type === 'function') ?? [];
    if (toolCalls.length > 0) {
      for (const toolCall of toolCalls) {
        const functionName = toolCall.function.name;
        const args = JSON.parse(toolCall.function.arguments);

        let result: string;

        try {
          switch (functionName) {
          case 'list_ui_elements':
            result = await context.listUIElements();
            options.aiPanel?.addUiMessage('Looking for the UI element', false);
            break;
          case 'get_available_viewers':
            result = context.getAvailableViewers();
            options.aiPanel?.addUiMessage('Getting available viewer types', false);
            break;
          case 'highlight_element':
            result = 'Element highlighted';
            const actRes = await context.highlightElement(args.elementId, args.description);
            options.aiPanel?.addUiMessage(`Highlighting: **${args.elementId}**`, false);
            if (actRes != null) {
              // if it is not null, then it is an error message
              result = actRes;
              options.aiPanel?.addUiMessage(`⚠️ ${actRes}`, false);
            }
            break;
          case 'describe_viewer':
            result = await context.describeViewer(args.viewerType);
            options.aiPanel?.addUiMessage(`Getting available properties for **${args.viewerType}** viewer`, false);
            break;
          case 'add_viewer': {
            const viewerProps = args.viewerProperties || args.properties || {};
            result = await context.addViewer(args.viewerType, viewerProps);
            options.aiPanel?.addUiMessage(`Added **${args.viewerType}** viewer ${viewerProps && Object.keys(viewerProps).length > 0 ? `with properties: \n\`\`\`json\n${JSON.stringify(viewerProps, null, 2)}\n\`\`\`\`` : ''}`, false);
            break;
          }
          case 'find_viewers_by_type':
            result = await context.findViewersByType(args.viewerType);
            options.aiPanel?.addUiMessage(`Searching for needed **${args.viewerType}** viewers`, false);
            break;
          case 'adjust_viewer': {
            const viewerProps = args.viewerProperties || args.properties || {};
            result = await context.adjustViewer(args.viewerId, viewerProps);
            options.aiPanel?.addUiMessage(`Setting properties to the viewer: \n\`\`\`json\n${JSON.stringify(viewerProps, null, 2)}\n\`\`\`\``, false);
            break;
          }
          case 'list_current_viewers':
            result = await context.listCurrentViewers();
            options.aiPanel?.addUiMessage('Listing currently open viewers', false);
            break;
          case 'reply_to_user':
            result = 'Message sent to user';
            options.aiPanel?.addUiMessage(args.message, false);
            break;
          default:
            result = `Unknown function: ${functionName}`;
          }
        } catch (error: any) {
          result = `Error: ${error.message}`;
          options.aiPanel?.addUiMessage(`⚠️ Error executing ${functionName}: ${error.message}`, false);
        }

        // Add tool result to messages
        const toolMessage: OpenAI.Chat.ChatCompletionToolMessageParam = {
          role: 'tool',
          content: result,
          tool_call_id: toolCall.id
        };
        messages.push(toolMessage);
        options.aiPanel?.addEngineMessage({...toolMessage});
      }
    } else {
      // No tool calls and has content - we're done
      if (assistantMessage.content)
        break;
    }

    // Check if we should continue
    if (response.choices[0].finish_reason === 'stop')
      break;

    // Check if we've hit max iterations
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

  async addViewer(viewerType: string, properties: Record<string, any>): Promise<string> {
    const newViewer = this.tableView.addViewer(viewerType, properties);
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

  async adjustViewer(viewerId: string, properties: any): Promise<string> {
    const foundViewer = Array.from(this.tableView.viewers).find((v) => this.getViewerIdentifier(v) === viewerId);
    if (!foundViewer)
      return `Viewer with ID ${viewerId} not found. Please add a new viewer instead.`;
    foundViewer.setOptions(properties);
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
  ] as const;

  async listUIElements(): Promise<string> {
    // TODO: Populate with actual UI elements from the table view
    const elements = this.highlightableElements.filter((el) => document.querySelector(el.selector) != null);

    let result = 'Available UI elements:\n\n';
    for (const elem of elements)
      result += `- **${elem.id}**: ${elem.description}\n`;

    return result;
  }

  async highlightElement(elementId: string, description: string): Promise<string | null> {
    const element = this.highlightableElements.find((el) => el.id.toLowerCase() === elementId.toLowerCase());
    if (!element)
      return `Element with ID ${elementId} not found or cannot be highlighted.`;
    const domElement = document.querySelector(element.selector) as HTMLElement;
    if (!domElement)
      return `Element with ID ${elementId} not found in the current UI.`;
    ui.hints.addHintIndicator(domElement, true, 10000);
    return null;
  }
}
