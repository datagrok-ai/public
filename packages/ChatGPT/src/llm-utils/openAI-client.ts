/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import OpenAI from 'openai';
import * as api from '../package-api';
import {LLMCredsManager} from './creds';
import {ChatCompletionParseParams} from 'openai/resources/chat/completions';
import {JsonSchema} from '../prompt-engine/interfaces';
export class OpenAIClient {
  public openai: OpenAI;
  private constructor(private vectorStoreId: string) {
    if (!grok.ai.openAiConfigured)
      throw new Error('OpenAI is not configured. Please set the API key in the AI settings under Admin section.');

    this.openai = new OpenAI({baseURL: `${grok.ai.openAiProxyUrl}/v1`, dangerouslyAllowBrowser: true,
      apiKey: 'unused',
      fetch: async (input, init) => {
        const url = typeof input === 'string' ? input : input.toString();
        return fetch(url, {
          ...init,
          headers: {
            ...init?.headers,
            Authorization: grok.ai.openAiProxyToken,
          },
        });
      },
    });
  }

  private static _instance: OpenAIClient | null = null;
  static getInstance(): OpenAIClient {
    OpenAIClient._instance ??= new OpenAIClient(LLMCredsManager.getVectorStoreId());
    return OpenAIClient._instance;
  }

  async getHelpAnswer(question: string): Promise<string> {
    if (!this.vectorStoreId)
      throw new Error('Vector Store ID is not set. Please set it in the package settings.');
    const storageIDs = [this.vectorStoreId];
    const response = await this.openai.responses.create({
      model: 'gpt-4o-mini',
      instructions: `
    You are a helpful assistant with access to documentation for the datagrok platform via file_search tool. Use the available tool to answer questions accurately.
    
    VERY IMPORTANT TO REMEMBER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, Always answer with UI solution first (if any) and then the code (if relevant).

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

  async generalPromptCached(model: string, systemPrompt: string, prompt: string, schema?: JsonSchema): Promise<string> {
    return await api.funcs.askAIGeneralCached(model, systemPrompt, prompt, schema);
  }

  /**
   * @deprecated - use generalPromptCached instead
   * @param model - model name
   * @param systemPrompt - system prompt
   * @param prompt - user prompt
   * @returns string response from OpenAI
   */
  async generalPrompt(
    model: string,
    systemPrompt: string,
    userPrompt: string,
    schema?: JsonSchema
  ): Promise<string> {
    const request: ChatCompletionParseParams = {
      model,
      temperature: model.startsWith('gpt-5') ? 1 : 0, // 5th models only work with temperature 1
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
      return JSON.stringify(msg?.parsed);
    return msg?.content ?? '';
  }
}

export async function askOpenAIHelp(question: string): Promise<string> {
  return await api.funcs.askDocumentationCached(question);
}
