/* eslint-disable max-len */
import OpenAI from 'openai';
import * as api from '../package-api';
import {LLMCredsManager} from './creds';
import { ChatCompletionParseParams } from 'openai/resources/chat/completions';
export class OpenAIHelpClient {
  private openai: OpenAI;
  private constructor(private apiKey: string, private vectorStoreId: string) {
    this.openai = new OpenAI({
      apiKey: this.apiKey, dangerouslyAllowBrowser: true
    });
  }

  private static _instance: OpenAIHelpClient | null = null;
  static getInstance(): OpenAIHelpClient {
    OpenAIHelpClient._instance ??= new OpenAIHelpClient(LLMCredsManager.getApiKey(), LLMCredsManager.getVectorStoreId());
    return OpenAIHelpClient._instance;
  }

  async getHelpAnswer(question: string): Promise<string> {
    const storageIDs = [this.vectorStoreId];
    const response = await this.openai.responses.create({
      model: 'gpt-4o',
      instructions: `
    You are a helpful assistant with access to documentation for the datagrok platform via file_search tool. Use the available tool to answer questions accurately.
    
    VERY IMPORTANT TO REMEMER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, Always answer with UI solution first (if any) and then the code (if relevant).

    Make sure the output is nicely formatted with markdown syntax where applicable.

    If the question is also related to coding, provide code examples in your answers.

    MAKE SURE TO ALWAYS USE TOOL to access the repository information! DO NOT MAKE UP ANSWERS BASED ON YOUR TRAINING DATA!

    Prioritize readme files from help directory, js file samples from apisamples directory, and any other relevant documentation files.
    `,
      input: `${question}.\n Remember that I am in Datagrok platform environment and the question is related to working in Datagrok platform.`,
      tools: [
        {
          type: 'file_search',
          vector_store_ids: storageIDs,
        },
      ],
    });
    return response.output_text;
  }

  async generalPromptCached(model: string, systemPrompt: string, prompt: string, schema?: { [key: string]: unknown }): Promise<string> {
    return await api.funcs.askAIGeneralCached(model, systemPrompt, prompt, schema);
  }

  /**
   * @deprecated - use generalPromptCached instead
   * @param model - model name
   * @param systemPrompt - system prompt
   * @param prompt - user prompt
   * @returns string response from OpenAI
   */
  async generalPrompt<T = any>(
    model: string,
    systemPrompt: string,
    userPrompt: string,
    schema?: Record<string, unknown>
  ): Promise<T | string> {
    const request: ChatCompletionParseParams = {
      model,
      temperature: 0,
      messages: [
        {role: 'system', content: systemPrompt},
        {role: 'user', content: userPrompt},
      ],
      ...(schema && {
        response_format: {
          type: 'json_schema',
          json_schema: {
            name: 'ResponseSchema',
            schema,
          },
        },
      }),
    };

    const response = await this.openai.chat.completions.parse(request);
    const msg = response.choices?.[0]?.message;

    if (schema)
      return msg?.parsed as T;
    return msg?.content ?? '';
  }
}

export async function askOpenAIHelp(question: string): Promise<string> {
  return await api.funcs.askDocumentationCached(question);
}
