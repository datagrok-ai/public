import * as DG from 'datagrok-api/dg';

const OPENAIApiKeyName = 'apiKey';
const OPENAIVectorStoreName = 'vectorStoreId';

export class LLMCredsManager {
  apiKey: string;
  vectorStoreId: string;
  private constructor(apiKey: string, vectorStoreId: string) {
    this.apiKey = apiKey;
    this.vectorStoreId = vectorStoreId;
  }

  private static _instance: LLMCredsManager | null = null;
  // Initialize the singleton instance
  // TODO: in future this will probably move to credentials manager
  static init(pckg: DG.Package): void {
    const apiKey = pckg.settings[OPENAIApiKeyName];
    const vectorStoreId = pckg.settings[OPENAIVectorStoreName];
    LLMCredsManager._instance = new LLMCredsManager(apiKey, vectorStoreId);
  }
  static getInstance(): LLMCredsManager {
    return LLMCredsManager._instance!;
  }

  static isAvailable(): boolean {
    const instance = LLMCredsManager.getInstance();
    return instance.apiKey.length > 0 && instance.vectorStoreId.length > 0;
  }

  static getApiKey(): string {
    return LLMCredsManager.getInstance().apiKey;
  }
  static getVectorStoreId(): string {
    return LLMCredsManager.getInstance().vectorStoreId;
  }
}
