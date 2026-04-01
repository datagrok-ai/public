export enum LLMApiNames {
  OpenAIChatCompletions = 'OpenAI Chat Completions',
  OpenAIResponses = 'OpenAI Responses',
  AntropicMessages = 'Antropic Messages',
  AmazonBedrock = 'Amazon Bedrock',
}

export type LLMApiName = `${LLMApiNames}`;

export type PackageSettings = {
  vectorStoreId?: string;
  defaultFastModel?: string;
  defaultDeepResearchModel?: string;
  defaultCodingModel?: string;
  APIName: LLMApiName;
  region?: string;
  apiVersion?: string;
}
