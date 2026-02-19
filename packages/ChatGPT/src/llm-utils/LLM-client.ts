/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import OpenAI from 'openai';
import * as api from '../package-api';
import {LLMCredsManager} from './creds';
import {JsonSchema, LLMApiNames} from '../prompt-engine/interfaces';
import {_package} from '../package';
// import {AIProvider, BuiltinToolSpec} from './AI-API-providers/types';
import * as osdk from '@ai-sdk/openai';
import * as asdk from '@ai-sdk/anthropic';
import {LanguageModelV3Message, LanguageModelV3} from '@ai-sdk/provider';
import {findLast} from '../utils';

export type ModelOption = 'Fast' | 'Deep Research' | 'Coding';
export const ModelType: {[type in ModelOption]: string} = {
  Fast: 'gpt-4o-mini',
  'Deep Research': 'gpt-5.2', // hell of a smart model but a bit expensive expensive
  Coding: 'gpt-5.1-codex-max',
  //'Deep Research': 'o4-mini', // good balance between speed, quality and $$$
};

export function initModelTypeNames() {
  if (_package.settings.defaultFastModel)
    ModelType.Fast = _package.settings.defaultFastModel;
  if (_package.settings.defaultDeepResearchModel)
    ModelType['Deep Research'] = _package.settings.defaultDeepResearchModel;
  if (_package.settings.defaultCodingModel)
    ModelType.Coding = _package.settings.defaultCodingModel;
  if (grok.ai.config.current?.provider == 'azure' && grok.ai.config.current.modelToDeployment) {
    const cfg: typeof grok.ai.config.current.modelToDeployment = DG.toJs(grok.ai.config.current.modelToDeployment);
    ModelType.Fast = cfg[ModelType.Fast] ?? ModelType.Fast;
    ModelType['Deep Research'] = cfg[ModelType['Deep Research']] ?? ModelType['Deep Research'];
    ModelType.Coding = cfg[ModelType.Coding] ?? ModelType.Coding;
  }
}

export class LLMClient {
  public openai: OpenAI;
  public aiModels!: Record<ModelOption, LanguageModelV3>;
  private constructor(private vectorStoreId: string) {
    if (!grok.ai.config.configured)
      throw new Error('AI is not configured. Please set the API key in the AI settings under Admin section.');
    const apiVersion = _package.settings.apiVersion ?? '';
    const versionChunk = apiVersion ? `/${apiVersion}` : '';
    const fullBaseUrl = `${grok.ai.config.proxyUrl}${versionChunk}`;
    // this is still needed for vector store access
    this.openai = new OpenAI({baseURL: fullBaseUrl, dangerouslyAllowBrowser: true,
      apiKey: 'unused',
      // ...(grok.ai.config.apiMode === 'legacy' ? {apiVersion: grok.ai.config!.apiVersion!} : {}),
      fetch: async (input, init) => {
        let url = typeof input === 'string' ? input : input.toString();
        if (grok.ai.config.current?.provider == 'azure' && grok.ai.config.current.apiMode === 'legacy') {
          if (grok.ai.config.current.apiVersion) {
            const prsed = new URL(url);
            prsed.searchParams.set('api-version', grok.ai.config.current.apiVersion);
            url = prsed.toString();
          }
          if (init?.body && typeof init.body === 'string') {
            const bodyObj = JSON.parse(init.body);
            const modelName = bodyObj['model'];
            if (!modelName)
              throw new Error('Model name is not specified in the request body.');
            const modifiedURL = url.substring(0, fullBaseUrl.length) + '/' + modelName + url.substring(fullBaseUrl.length);
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
        headers.set('Authorization', grok.ai.config.proxyToken);
        prms.headers = headers;
        return fetch(url, prms);
      },
    });

    // const constr = _package.settings.APIName === 'openai chat completions' ? osdk.openai.chat : osdk.openai.responses;
    // this.aiModels.Fast = constr(ModelType.Fast);
    const providerConstructor = _package.settings.APIName === LLMApiNames.AntropicMessages ? asdk.createAnthropic : osdk.createOpenAI;
    const provider = providerConstructor({
      baseURL: fullBaseUrl,
      apiKey: grok.ai.config.proxyToken!,
      fetch: async (input, init) => {
        let url = typeof input === 'string' ? input : input.toString();
        const currentConfig = grok.ai.config.current;
        if ('apiMode' in currentConfig && currentConfig.apiMode === 'legacy' && currentConfig.apiVersion) {
          const prsed = new URL(url);
          prsed.searchParams.set('api-version', currentConfig.apiVersion);
          url = prsed.toString();
          if (init?.body && typeof init.body === 'string') {
            const bodyObj = JSON.parse(init.body);
            const modelName = bodyObj['model'];
            if (!modelName)
              throw new Error('Model name is not specified in the request body.');
            const modifiedURL = url.substring(0, fullBaseUrl.length) + '/' + modelName + url.substring(fullBaseUrl.length);
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
        headers.set('Authorization', grok.ai.config.proxyToken);
        prms.headers = headers;
        return fetch(url, prms);
      },
    });

    const constr = _package.settings.APIName === LLMApiNames.OpenAIChatCompletions ? provider.chat :
      (_package.settings.APIName === LLMApiNames.OpenAIResponses ? (provider as osdk.OpenAIProvider).responses :
        (provider as asdk.AnthropicProvider).messages);
    // @ts-ignore
    this.aiModels = {};
    for (const type of Object.keys(ModelType) as ModelOption[])
      this.aiModels[type] = constr(ModelType[type]);
    // this.aiModels.Fast = constr(ModelType.Fast);
    // this.aiModels['Deep Research'] = constr(ModelType['Deep Research']);
    // this.aiModels.Coding = constr(ModelType.Coding);
    // this.aiModels.Fast.doGenerate({prompt: [{role: 'user', content: [{type: 'text', text: 'call the test function with param1 equal to hello'}]}],
    //   tools: [{
    //     type: 'function', name: 'test', description: 'Make sure to call the test function', inputSchema: {type: 'object', properties: {param1: {type: 'string', description: 'a test param'}}, required: ['param1']}
    //   }]})
    //   .then((res) => { console.log('OpenAI client initialized', res); });
  }

  private static _instance: LLMClient | null = null;
  static getInstance(): LLMClient {
    LLMClient._instance ??= new LLMClient(LLMCredsManager.getVectorStoreId());
    return LLMClient._instance;
  }

  async getHelpAnswer(question: string): Promise<string> {
    if (!this.vectorStoreId)
      throw new Error('Vector Store ID is not set. Please set it in the package settings.');
    const storageIDs = [this.vectorStoreId];
    // const provider = AIProvider.getProvider();
    const systemMsg: LanguageModelV3Message = {
      role: 'system',
      content: `
    You are a helpful assistant with access to documentation for the datagrok platform via file_search tool. Use the available tool to answer questions accurately.
    
    VERY IMPORTANT TO REMEMBER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, Always answer with UI solution first (if any) and then the code (if relevant).

    Make sure the output is nicely formatted with markdown syntax where applicable.

    If the question is also related to coding, provide code examples in your answers.

    MAKE SURE TO ALWAYS USE TOOL to access the repository information! DO NOT MAKE UP ANSWERS BASED ON YOUR TRAINING DATA!

    Prioritize readme files from help directory, js file samples from apisamples directory, and any other relevant documentation files.
    `
    };
    const userMsg: LanguageModelV3Message = {
      role: 'user',
      content: [{type: 'text', text: `${question}.\n Remember that I am in Datagrok platform environment and the question is related to working in Datagrok platform.`}]
    };
    // const tools: BuiltinToolSpec[] = [{type: 'file_search', vector_store_ids: storageIDs}];
    // const response = await provider.create({
    //   model: ModelType.Fast,
    //   messages: [systemMsg, userMsg],
    //   tools,
    // });

    const response = await this.aiModels.Fast.doGenerate({
      prompt: [systemMsg, userMsg],
      tools: [{type: 'provider', id: 'openai.file_search', args: {vectorStoreIds: storageIDs}, name: 'file_search'}],
    });

    const output = response.content.length > 0 ? findLast(response.content, (c) => c.type === 'text') : undefined;

    return output?.text ?? '';
  }

  async generalPromptCached(model: string, systemPrompt: string, prompt: string, schema?: JsonSchema): Promise<string> {
    return await api.funcs.askAIGeneralCached(model, systemPrompt, prompt, schema);
  }

  public createSystemMessage(content: string): LanguageModelV3Message {
    return {
      role: 'system',
      content: content,
    };
  }

  public createToolOutputMessage(callId: string, toolName: string, output: string): LanguageModelV3Message {
    return {
      role: 'tool',
      content: [{type: 'tool-result', toolCallId: callId, toolName: toolName, output: {type: 'text', value: output}}],
    };
  }

  public createUserMessage(content: string): LanguageModelV3Message {
    return {
      role: 'user',
      content: [{type: 'text', text: content}],
    };
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
    const systemMsg = this.createSystemMessage(systemPrompt);
    const userMsg = this.createUserMessage(userPrompt);

    const modelProvider = Object.values(this.aiModels).find((m) => m.modelId === model) ?? this.aiModels.Fast;
    const res = await modelProvider.doGenerate({
      prompt: [systemMsg, userMsg],
      responseFormat: {
        ...(!schema ? {type: 'text'} : {type: 'json', name: 'ResponseSchema', schema: schema!}),
      }
    });

    const output = res.content.length > 0 ? findLast(res.content, (c) => c.type === 'text') : undefined;
    const text = output?.text ?? '';
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
