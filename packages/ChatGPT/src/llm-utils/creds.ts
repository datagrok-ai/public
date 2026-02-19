import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

const OPENAIVectorStoreName = 'vectorStoreId';

export class LLMCredsManager {
  vectorStoreId: string;
  private constructor(vectorStoreId: string) {
    this.vectorStoreId = vectorStoreId;
  }

  private static _instance: LLMCredsManager | null = null;
  // Initialize the singleton instance
  // TODO: in future this will probably move to credentials manager
  static init(pckg: DG.Package): void {
    const hasAIKey = grok.ai.config.configured;
    const vectorStoreId = pckg.settings[OPENAIVectorStoreName];
    LLMCredsManager._instance =
      new LLMCredsManager(hasAIKey ? vectorStoreId ?? '' : '');
  }
  static getInstance(): LLMCredsManager {
    return LLMCredsManager._instance!;
  }
  static getVectorStoreId(): string {
    const res = LLMCredsManager.getInstance().vectorStoreId;
    if (!res || !res.trim())
      return '';
    return res;
  }
}
