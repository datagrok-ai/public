// /* eslint-disable max-len */

// export type AIProviderEndpoint =
//     | 'openai.responses'
//     | 'openai.chat.completions';

// export type AIProviderToolCall = {
//     id: string;
//     type: 'function';
//     name: string;
//     arguments: string;
// };

// export type AIProviderUsage = {
//     inputTokens?: number;
//     outputTokens?: number;
//     totalTokens?: number;
// };

// export type AIProviderReasoning = {
//     effort?: 'low' | 'medium' | 'high';
// };

// export type EmptyObject = Record<string, never>;

// export type JsonSchema = {
//     type: 'string' | 'array' | 'object' | 'number' | 'integer' | 'boolean';
//     description?: string;
//     items?: JsonSchema;
//     properties?: Record<string, JsonSchema>;
//     required?: string[];
//     additionalProperties?: boolean | JsonSchema;
//     patternProperties?: Record<string, JsonSchema>;
//     default?: string | number | boolean | null;
// };

// export type FunctionToolSpec = {
//     type: 'function';
//     name: string;
//     description: string;
//     parameters: JsonSchema;
//     strict: true;
// };

// export type BuiltinToolSpec = {
//     type: 'file_search';
//     vector_store_ids: string[];
// };

// export type ToolSpec = FunctionToolSpec | BuiltinToolSpec;

// export type AIProviderRequest<
//     TModel,
//     TInputMessage,
//     TTool,
//     TToolChoice = never,
//     TResponseFormat = never,
//     TRequestExtra = EmptyObject
// > = {
//     model: TModel;
//     messages: TInputMessage[];
//     tools?: TTool[];
//     toolChoice?: TToolChoice;
//     temperature?: number;
//     maxOutputTokens?: number;
//     reasoning?: AIProviderReasoning;
//     responseFormat?: TResponseFormat;
//     extra?: TRequestExtra;
// };

// export type AIProviderResult<TMessage, TRawResponse> = {
//     raw: TRawResponse;
//     outputText: string;
//     outputMessages: TMessage[];
//     toolCalls: AIProviderToolCall[];
//     usage?: AIProviderUsage;
// };

// export interface AIAPIProvider<
//     TModel,
//     TMessage,
//     TTool,
//     TRawResponse,
//     TToolChoice = never,
//     TResponseFormat = never,
//     TRequestExtra = EmptyObject
// > {
//     readonly endpoint: AIProviderEndpoint;
//     create(request: AIProviderRequest<TModel, TMessage, TTool, TToolChoice, TResponseFormat, TRequestExtra>): Promise<AIProviderResult<TMessage, TRawResponse>>;
//     createSystemMessage(content: string): TMessage;
//     createUserMessage(content: string): TMessage;
//     createToolOutputMessage(callId: string, output: string): TMessage;
// }

// export class AIProvider {
//   private static provider: AIAPIProvider<any, any, any, any, any, any, any> | null = null;
//   public static setProvider(provider: AIAPIProvider<any, any, any, any, any, any, any>): void {
//     AIProvider.provider = provider;
//   }
//   public static getProvider(): AIAPIProvider<any, any, any, any, any, any, any> {
//     if (!AIProvider.provider)
//       throw new Error('AI Provider is not set');

//     return AIProvider.provider;
//   }
// }
