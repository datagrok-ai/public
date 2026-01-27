/* eslint-disable max-len */
import type OpenAI from 'openai';
import type {ChatModel} from 'openai/resources/shared';
import type {
  AIAPIProvider,
  AIProviderRequest,
  AIProviderResult,
  AIProviderToolCall,
  ToolSpec,
} from './types';
import type {MessageType} from '../panel';

function isObjectRecord(value: object | null | undefined): value is Record<string, string | number | boolean | null | object> {
  return !!value && typeof value === 'object';
}

type ResponseCreateBody = OpenAI.Responses.ResponseCreateParamsNonStreaming;
export type OpenAIResponsesTool = Exclude<ResponseCreateBody['tools'], undefined>[number];
export type OpenAIResponsesInput = Exclude<ResponseCreateBody['input'], string | undefined>;
export type OpenAIResponsesOutput = OpenAI.Responses.ResponseOutputItem;

export class OpenAIResponsesProvider
implements AIAPIProvider<
  ChatModel,
  MessageType,
  ToolSpec,
  OpenAI.Responses.Response,
  ResponseCreateBody['tool_choice'],
  OpenAI.Responses.ResponseTextConfig | OpenAI.Responses.ResponseTextConfig['format'],
  Partial<ResponseCreateBody>
> {
  readonly endpoint = 'openai.responses' as const;

  constructor(private readonly openai: OpenAI) {}

  createSystemMessage(content: string): MessageType {
    return {role: 'system', content} as OpenAI.Responses.ResponseInputItem;
  }

  createUserMessage(content: string): MessageType {
    return {role: 'user', content} as OpenAI.Responses.ResponseInputItem;
  }

  createToolOutputMessage(callId: string, output: string): MessageType {
    return {type: 'function_call_output', call_id: callId, output} as OpenAI.Responses.ResponseInputItem;
  }

  async create(
    request: AIProviderRequest<
      ChatModel,
      MessageType,
      ToolSpec,
      ResponseCreateBody['tool_choice'],
      OpenAI.Responses.ResponseTextConfig | OpenAI.Responses.ResponseTextConfig['format'],
      Partial<ResponseCreateBody>
    >
  ): Promise<AIProviderResult<MessageType, OpenAI.Responses.Response>> {
    const textConfig = request.responseFormat as
      | OpenAI.Responses.ResponseTextConfig
      | OpenAI.Responses.ResponseTextConfig['format']
      | undefined;

    const body: ResponseCreateBody = {
      model: request.model,
      input: request.messages as OpenAIResponsesInput,
      ...(request.tools ? {tools: request.tools} : {}),
      ...(request.toolChoice ? {tool_choice: request.toolChoice as ResponseCreateBody['tool_choice']} : {}),
      ...(typeof request.temperature === 'number' ? {temperature: request.temperature} : {}),
      ...(typeof request.maxOutputTokens === 'number' ? {max_output_tokens: request.maxOutputTokens} : {}),
      ...(request.reasoning ? {reasoning: request.reasoning as ResponseCreateBody['reasoning']} : {}),
      ...(request.extra ?? {}),
      stream: false,
    };

    if (textConfig) {
      if (isObjectRecord(textConfig) && ('format' in textConfig || 'verbosity' in textConfig))
        body.text = textConfig as OpenAI.Responses.ResponseTextConfig;
      else
        body.text = {format: textConfig as OpenAI.Responses.ResponseTextConfig['format']};
    }

    const response = await this.openai.responses.create(body);

    const outputItems = response.output ?? [];
    const toolCalls: AIProviderToolCall[] = [];

    for (const item of outputItems) {
      if (item.type === 'function_call') {
        const id = typeof item.call_id === 'string' ? item.call_id : '';
        const name = typeof item.name === 'string' ? item.name : '';
        const args = typeof item.arguments === 'string' ? item.arguments : '{}';
        toolCalls.push({id, type: 'function', name, arguments: args});
      }
    }

    return {
      raw: response,
      outputText: (response.output_text ?? '').trim(),
      outputMessages: outputItems as MessageType[],
      toolCalls,
      usage: response.usage ? {
        inputTokens: response.usage.input_tokens ?? undefined,
        outputTokens: response.usage.output_tokens ?? undefined,
        totalTokens: response.usage.total_tokens ?? undefined,
      } : undefined,
    };
  }
}
