// /* eslint-disable max-len */
// import type OpenAI from 'openai';
// import type {ChatModel} from 'openai/resources/shared';
// import type {
//   AIAPIProvider,
//   AIProviderRequest,
//   AIProviderResult,
//   AIProviderToolCall,
//   ToolSpec,
// } from './types';
// import type {MessageType} from '../panel';

// type ChatCreateBody = OpenAI.Chat.ChatCompletionCreateParamsNonStreaming;
// export type OpenAIChatTool = Exclude<ChatCreateBody['tools'], undefined>[number];
// export type OpenAIChatMessage = OpenAI.Chat.ChatCompletionMessageParam;

// export class OpenAIChatCompletionsProvider
// implements AIAPIProvider<
//   ChatModel,
//   MessageType,
//   ToolSpec,
//   OpenAI.Chat.ChatCompletion,
//   ChatCreateBody['tool_choice'],
//   ChatCreateBody['response_format'],
//   Partial<ChatCreateBody>
// > {
//   readonly endpoint = 'openai.chat.completions' as const;

//   constructor(private readonly openai: OpenAI) {}

//   createSystemMessage(content: string): MessageType {
//     return {role: 'system', content} as OpenAI.Chat.ChatCompletionMessageParam;
//   }

//   createUserMessage(content: string): MessageType {
//     return {role: 'user', content} as OpenAI.Chat.ChatCompletionMessageParam;
//   }

//   createToolOutputMessage(callId: string, output: string): MessageType {
//     return {role: 'tool', tool_call_id: callId, content: output} as OpenAI.Chat.ChatCompletionToolMessageParam;
//   }

//   async create(
//     request: AIProviderRequest<
//       ChatModel,
//       MessageType,
//       ToolSpec,
//       ChatCreateBody['tool_choice'],
//       ChatCreateBody['response_format'],
//       Partial<ChatCreateBody>
//     >
//   ): Promise<AIProviderResult<MessageType, OpenAI.Chat.ChatCompletion>> {
//     const tools = request.tools?.map((tool) => {
//       if (tool.type !== 'function')
//         throw new Error(`Chat Completions does not support tool type: ${tool.type}`);
//       return {
//         type: 'function',
//         function: {
//           name: tool.name,
//           description: tool.description,
//           parameters: tool.parameters,
//           strict: tool.strict,
//         },
//       } as OpenAIChatTool;
//     });

//     const body: ChatCreateBody = {
//       model: request.model,
//       messages: request.messages as OpenAIChatMessage[],
//       ...(tools ? {tools} : {}),
//       ...(request.toolChoice ? {tool_choice: request.toolChoice as ChatCreateBody['tool_choice']} : {}),
//       ...(typeof request.temperature === 'number' ? {temperature: request.temperature} : {}),
//       ...(typeof request.maxOutputTokens === 'number' ? {max_completion_tokens: request.maxOutputTokens} : {}),
//       ...(request.reasoning?.effort ? {reasoning_effort: request.reasoning.effort} : {}),
//       ...(request.responseFormat ? {response_format: request.responseFormat as ChatCreateBody['response_format']} : {}),
//       ...(request.extra ?? {}),
//       stream: false,
//     };

//     const response = await this.openai.chat.completions.create(body);

//     const outputMessages = response.choices.map((choice: OpenAI.Chat.ChatCompletion.Choice) => choice.message);
//     const toolCalls: AIProviderToolCall[] = [];

//     for (const message of outputMessages) {
//       if (!message || typeof message !== 'object')
//         continue;
//       const calls = (message as OpenAI.Chat.ChatCompletionMessage).tool_calls;
//       if (Array.isArray(calls)) {
//         for (const call of calls) {
//           if (call.type === 'function') {
//             const id = call.id ?? '';
//             const name = call.function?.name ?? '';
//             const args = call.function?.arguments ?? '{}';
//             toolCalls.push({id, type: 'function', name, arguments: args});
//           }
//         }
//       }
//     }

//     const firstContent = outputMessages.find((m: OpenAI.Chat.ChatCompletionMessage) => typeof m?.content === 'string')?.content ?? '';

//     return {
//       raw: response,
//       outputText: firstContent.trim(),
//       outputMessages: outputMessages as MessageType[],
//       toolCalls,
//       usage: response.usage ? {
//         inputTokens: response.usage.prompt_tokens ?? undefined,
//         outputTokens: response.usage.completion_tokens ?? undefined,
//         totalTokens: response.usage.total_tokens ?? undefined,
//       } : undefined,
//     };
//   }
// }
