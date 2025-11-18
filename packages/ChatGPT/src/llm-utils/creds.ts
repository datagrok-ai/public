import * as DG from 'datagrok-api/dg';

export const OPENAI_HELP_APIKEY_STORAGE_NAME = 'openAIHelpApiKey';
export const OPENAI_HELP_STORAGE_ID = 'openAIHelpVectorStoreId';

export class OpenAIHelpClientCred {
  apiKey: string;
  vectorStoreId: string;
  private constructor(apiKey: string, vectorStoreId: string) {
    this.apiKey = apiKey;
    this.vectorStoreId = vectorStoreId;
  }

  private static _instance: OpenAIHelpClientCred | null = null;
  static async init(pckg: DG.Package): Promise<void> {
    const creds = (await pckg.getCredentials())?.parameters;
    if (creds) {
      const apiKey = creds[OPENAI_HELP_APIKEY_STORAGE_NAME] as string;
      const vectorStoreId = creds[OPENAI_HELP_STORAGE_ID] as string;
      if (apiKey && vectorStoreId)
        OpenAIHelpClientCred._instance ??= new OpenAIHelpClientCred(apiKey, vectorStoreId);
      else
        console.warn('OpenAI Help Client credentials are not set properly.');
    }
  }
  static getInstance(): OpenAIHelpClientCred | null {
    return OpenAIHelpClientCred._instance;
  }
}
