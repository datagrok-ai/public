/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import OpenAI from 'openai';
import * as api from '../package-api';
import {LLMCredsManager} from './creds';
import {JsonSchema} from '../prompt-engine/interfaces';
import {ChatModel} from 'openai/resources/shared';
import {_package} from '../package';
import {AIProvider, BuiltinToolSpec} from './AI-API-providers/types';
import {AzureOpenAI} from 'openai';

export type ModelOption = 'Fast' | 'Deep Research' | 'Coding';
export const ModelType: {[type in ModelOption]: ChatModel} = {
  Fast: 'gpt-4o-mini',
  'Deep Research': 'gpt-5.2' as ChatModel, // hell of a smart model but a bit expensive expensive
  Coding: 'gpt-5.1-codex-max' as ChatModel,
  //'Deep Research': 'o4-mini', // good balance between speed, quality and $$$
};

export function initModelTypeNames() {
  if (_package.settings.defaultFastModel)
    ModelType.Fast = _package.settings.defaultFastModel;
  if (_package.settings.defaultDeepResearchModel)
    ModelType['Deep Research'] = _package.settings.defaultDeepResearchModel as ChatModel;
  if (_package.settings.defaultCodingModel)
    ModelType.Coding = _package.settings.defaultCodingModel as ChatModel;
  if (grok.ai.config?.modelToDeployment) {
    const cfg: typeof grok.ai.config.modelToDeployment = DG.toJs(grok.ai.config.modelToDeployment);
    ModelType.Fast = cfg[ModelType.Fast] as ChatModel ?? ModelType.Fast;
    ModelType['Deep Research'] = cfg[ModelType['Deep Research']] as ChatModel ?? ModelType['Deep Research'];
    ModelType.Coding = cfg[ModelType.Coding] as ChatModel ?? ModelType.Coding;
  }
}

export class OpenAIClient {
  public openai: OpenAI;
  private constructor(private vectorStoreId: string) {
    if (!grok.ai.openAiConfigured || !grok.ai.openAiProxyToken)
      throw new Error('OpenAI is not configured. Please set the API key in the AI settings under Admin section.');
    const apiVersion = _package.settings.apiVersion ?? '';
    const versionChunk = apiVersion ? `/${apiVersion}` : '';
    this.openai = new OpenAI({baseURL: `${grok.ai.openAiProxyUrl}${versionChunk}`, dangerouslyAllowBrowser: true,
      apiKey: 'unused',
      // ...(grok.ai.config.apiMode === 'legacy' ? {apiVersion: grok.ai.config!.apiVersion!} : {}),
      fetch: async (input, init) => {
        let url = typeof input === 'string' ? input : input.toString();
        if (grok.ai.config?.apiMode === 'legacy' && grok.ai.config?.apiVersion) {
          const prsed = new URL(url);
          prsed.searchParams.set('api-version', grok.ai.config?.apiVersion);
          url = prsed.toString();
          if (init?.body && typeof init.body === 'string') {
            const bodyObj = JSON.parse(init.body);
            const modelName = bodyObj['model'];
            if (!modelName)
              throw new Error('Model name is not specified in the request body.');
            const modifiedURL = url.substring(0, grok.ai.openAiProxyUrl.length) + '/' + modelName + url.substring(grok.ai.openAiProxyUrl.length);
            url = modifiedURL;
          }
        }
        const prms = {
          ...init,
          // headers: {
          //   ...in
          //   Authorization: grok.ai.openAiProxyToken,
          // },
        };
        const headers = new Headers(prms.headers);
        headers.set('Authorization', grok.ai.openAiProxyToken);
        prms.headers = headers;
        return fetch(url, prms);
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
    const provider = AIProvider.getProvider();

    const systemMsg = provider.createSystemMessage(`
    You are a helpful assistant with access to documentation for the datagrok platform via file_search tool. Use the available tool to answer questions accurately.
    
    VERY IMPORTANT TO REMEMBER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, Always answer with UI solution first (if any) and then the code (if relevant).

    Make sure the output is nicely formatted with markdown syntax where applicable.

    If the question is also related to coding, provide code examples in your answers.

    MAKE SURE TO ALWAYS USE TOOL to access the repository information! DO NOT MAKE UP ANSWERS BASED ON YOUR TRAINING DATA!

    Prioritize readme files from help directory, js file samples from apisamples directory, and any other relevant documentation files.
    `);
    const userMsg = provider.createUserMessage(`${question}.\n Remember that I am in Datagrok platform environment and the question is related to working in Datagrok platform.`);
    const tools: BuiltinToolSpec[] = [{type: 'file_search', vector_store_ids: storageIDs}];
    const response = await provider.create({
      model: ModelType.Fast,
      messages: [systemMsg, userMsg],
      tools,
    });
    return response.outputText;
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
    const provider = AIProvider.getProvider();
    const systemMsg = provider.createSystemMessage(systemPrompt);
    const userMsg = provider.createUserMessage(userPrompt);
    const responseFormat = schema ?
      (provider.endpoint === 'openai.responses' ?
        {format: {type: 'json_schema', name: 'ResponseSchema', schema, strict: true}} :
        {type: 'json_schema', json_schema: {name: 'ResponseSchema', schema, strict: true}}) :
      undefined;

    const response = await provider.create({
      model,
      messages: [systemMsg, userMsg],
      ...(responseFormat ? {responseFormat} : {}),
    });

    const text = response.outputText ?? '';
    if (!schema)
      return text;

    try {
      return JSON.stringify(JSON.parse(text));
    } catch {
      return text;
    }
  }
}

export async function askOpenAIHelp(question: string): Promise<string> {
  return await api.funcs.askDocumentationCached(question);
}
