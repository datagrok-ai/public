
// /** @deprecated NOT IN USE */
// class DeepWikiOpenAIClient {
//   private mcpClient: Client;
//   private openai: OpenAI;
//   private transport: StreamableHTTPClientTransport | null = null;
//   private repository: string;
//   private model: string;
//   private conversationHistory: OpenAI.Chat.ChatCompletionMessageParam[] = [];
//   private cachedTools: OpenAI.Chat.ChatCompletionTool[] | null = null;

//   private static _instance: DeepWikiOpenAIClient | null = null;
//   static getInstance(config: DeepWikiClientConfig): DeepWikiOpenAIClient {
//     DeepWikiOpenAIClient._instance ??= new DeepWikiOpenAIClient(config);
//     return DeepWikiOpenAIClient._instance;
//   }
//   private constructor(config: DeepWikiClientConfig) {
//     this.repository = config.repository;
//     this.model = config.openaiModel || 'gpt-4o-mini';

//     const apiKey = config.openaiApiKey!;
//     if (!apiKey)
//       throw new Error('OpenAI API key is required');


//     this.openai = new OpenAI({apiKey, dangerouslyAllowBrowser: true});

//     this.mcpClient = new Client(
//       {
//         name: 'deepwiki-openai-client',
//         version: '1.0.0',
//       },
//       {
//         capabilities: {},
//       }
//     );
//   }

//   async connect(): Promise<void> {
//     if (this.transport) {
//       console.log('Already connected to MCP server');
//       return;
//     }
//     console.log('Connecting to DeepWiki MCP server...');

//     // Connect to DeepWiki's free mcp endpoint
//     this.transport = new StreamableHTTPClientTransport(
//       new URL('https://mcp.deepwiki.com/mcp')
//     );

//     await this.mcpClient.connect(this.transport);

//     // List available tools
//     const toolsResponse = await this.mcpClient.listTools();
//     console.log(
//       '\n✓ Connected! Available tools:',
//       toolsResponse.tools.map((tool) => tool.name).join(', ')
//     );
//   }

//   async disconnect(): Promise<void> {
//     if (this.transport) {
//       await this.transport.close();
//       console.log('\n✓ Disconnected from DeepWiki MCP server');
//     }
//   }

//   private _askQuestionTool: OpenAI.Chat.Completions.ChatCompletionTool | null = null;
//   async askUsingAskTool(question: string): Promise<string> {
//     if (!this.transport)
//       throw new Error('Not connected to MCP server. Call connect() first.');

//     if (!this._askQuestionTool) {
//       const tools = await this.mcpClient.listTools();
//       const askTool = tools.tools.find((tool) => tool.name.includes('question'));
//       if (!askTool)
//         throw new Error('ask_question tool not found on MCP server.');

//       this._askQuestionTool = {
//         type: 'function',
//         function: {
//           name: askTool.name,
//           description: askTool.description || '',
//           parameters: askTool.inputSchema as Record<string, unknown>,
//         },
//       };
//     }

//     const q = `You are a helpful assistant with access to documentation for the ${this.repository} GitHub repository.
      
//       VERY IMPORTANT TO REMEMER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, try to answer with UI in mind and also the code.

//       Make sure the output is nicely formatted with markdown syntax where applicable.

//       Do not repeat the question in your answer, and do not include the deepwiki links!!!.

//       user: ${question}
//       `;
//     // just call the ask tool directly
//     const toolArgs = {
//       repoName: this.repository,
//       repository: this.repository,
//       question: q,
//     };

//     const result = await this.mcpClient.callTool({
//       name: this._askQuestionTool.function.name,
//       arguments: toolArgs,
//     });

//     const contentArray = Array.isArray(result.content) ? result.content : [];
//     const toolResult = contentArray
//       .filter((item: any) => item.type === 'text')
//       .map((item: any) => item.text || '')
//       .join('\n');
//     if (!toolResult || toolResult.length === 0)
//       return 'No response generated.';

//     // remove deepwiki links from the answer
//     const lines = toolResult.split('\n').filter((line) => !line.includes('deepwiki.com')).join('\n');
//     return lines;
//   }

//   async askQuestion(question: string): Promise<string> {
//     if (!this.transport)
//       throw new Error('Not connected to MCP server. Call connect() first.');


//     console.log(`\n🤔 Question: ${question}`);

//     // Get or use cached tools from MCP server
//     if (!this.cachedTools) {
//       const toolsResponse = await this.mcpClient.listTools();
//       this.cachedTools = toolsResponse.tools
//       //.filter((a) => !a.name.includes('question'))
//         .map(
//           (tool) => ({
//             type: 'function',
//             function: {
//               name: tool.name,
//               description: tool.description || '',
//               parameters: tool.inputSchema as Record<string, unknown>,
//             },
//           })
//         );
//     }

//     // System message (only sent once per conversation in history)
//     const systemMessage: OpenAI.Chat.ChatCompletionMessageParam = {
//       role: 'system',
//       content: `You are a helpful assistant with access to documentation for the ${this.repository} GitHub repository via DeepWiki MCP tools. Use the available tools to answer questions accurately.
      
//       VERY IMPORTANT TO REMEMER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, try to answer with UI in mind and also the code.

//       Make sure the output is nicely formatted with markdown syntax where applicable.

//       If you use the ask_question tool and it returns comprehensive answer, make sure not to warp it and return as is.

//       MAKE SURE TO ALWAYS USE TOOLS to access the repository information! DO NOT MAKE UP ANSWERS BASED ON YOUR TRAINING DATA!
//       `,
//     };

//     // Add user question to conversation history
//     this.conversationHistory.push({
//       role: 'user',
//       content: question,
//     });

//     // Initial call to OpenAI
//     let response = await this.openai.chat.completions.create({
//       model: this.model,
//       messages: [systemMessage, ...this.conversationHistory],
//       tools: this.cachedTools,
//       tool_choice: 'auto',
//     });

//     let assistantMessage = response.choices[0].message;
//     this.conversationHistory.push(assistantMessage);

//     // Handle tool calls
//     while (assistantMessage.tool_calls && assistantMessage.tool_calls.length > 0) {
//       console.log(`\n🔧 Calling ${assistantMessage.tool_calls.length} tool(s)...`);

//       for (const toolCall of assistantMessage.tool_calls) {
//         const toolName = toolCall.function.name;
//         const toolArgs = JSON.parse(toolCall.function.arguments);

//         // Add repository to the arguments
//         toolArgs.repository = this.repository;

//         console.log(`   → ${toolName}(${JSON.stringify(toolArgs)})`);

//         // Call MCP tool
//         const result = await this.mcpClient.callTool({
//           name: toolName,
//           arguments: toolArgs,
//         });

//         // Extract text content from result with proper type handling
//         const contentArray = Array.isArray(result.content) ? result.content : [];
//         const toolResult = contentArray
//           .filter((item: any) => item.type === 'text')
//           .map((item: any) => item.text || '')
//           .join('\n');

//         // Add tool result to conversation
//         this.conversationHistory.push({
//           role: 'tool',
//           tool_call_id: toolCall.id,
//           content: toolResult,
//         });
//       }

//       // Get next response from OpenAI
//       response = await this.openai.chat.completions.create({
//         model: this.model,
//         messages: [systemMessage, ...this.conversationHistory],
//         tools: this.cachedTools,
//       });

//       assistantMessage = response.choices[0].message;
//       this.conversationHistory.push(assistantMessage);
//     }

//     const finalAnswer = assistantMessage.content || 'No response generated.';
//     console.log(`\n💡 Answer: ${finalAnswer}`);

//     return finalAnswer;
//   }

//   clearHistory(): void {
//     this.conversationHistory = [];
//     console.log('\n✓ Conversation history cleared');
//   }
// }