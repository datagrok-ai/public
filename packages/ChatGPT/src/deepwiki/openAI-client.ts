/* eslint-disable max-len */
import OpenAI from 'openai';

class OpenAIHelpClient {
  private openai: OpenAI;
  private constructor(private apiKey: string, private vectorStoreId: string) {
    this.openai = new OpenAI({
      apiKey: this.apiKey, dangerouslyAllowBrowser: true,
    });
  }

  private static _instance: OpenAIHelpClient | null = null;
  static getInstance(apiKey: string, vectorStoreId: string): OpenAIHelpClient {
    OpenAIHelpClient._instance ??= new OpenAIHelpClient(apiKey, vectorStoreId);
    return OpenAIHelpClient._instance;
  }

  async getHelpAnswer(question: string): Promise<string> {
    const storageIDs = [this.vectorStoreId];
    const response = await this.openai.responses.create({
      model: 'gpt-4o-mini',
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
}

export async function askOpenAIHelp(question: string, apiKey: string, vectorStoreId: string): Promise<string> {
  if (!apiKey || !vectorStoreId)
    return 'API Key or Vector Store ID is not set. Please configure them in the ChatGPT package settings.';

  const client = OpenAIHelpClient.getInstance(apiKey, vectorStoreId);
  return await client.getHelpAnswer(question);
}
