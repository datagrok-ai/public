/* eslint-disable max-len */
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {SSEClientTransport} from '@modelcontextprotocol/sdk/client/sse';

interface DeepWikiClientConfig {
  repository: string;
  openaiApiKey: string;
  openaiModel?: string;
}

export const helpAIName = 'WIKITOR GROKOWSKY';

class DeepWikiMCPClient {
  private mcpClient: Client;
  private transport: SSEClientTransport | null = null;
  private readonly askTool = {
    name: 'ask_question',
    description: 'Ask a question about the repository and get an answer based on its documentation.',
  };
  private static _instance: DeepWikiMCPClient | null = null;
  static getInstance(repo?: string): DeepWikiMCPClient {
    DeepWikiMCPClient._instance ??= new DeepWikiMCPClient(repo);
    return DeepWikiMCPClient._instance;
  }
  private constructor(private repoName: string = 'datagrok-ai/public') {
    this.mcpClient = new Client(
      {
        name: 'deepwiki-GROKOWSKY-client',
        version: '1.0.0',
      },
      {
        capabilities: {},
      }
    );
  }

  async connect(): Promise<void> {
    if (this.transport) {
      console.log('Already connected to MCP server');
      return;
    }
    console.log('Connecting to DeepWiki MCP server...');

    // Connect to DeepWiki's free mcp endpoint
    this.transport = new SSEClientTransport(
      new URL('https://mcp.deepwiki.com/mcp')
    );

    await this.mcpClient.connect(this.transport);
  }

  async disconnect(): Promise<void> {
    if (this.transport) {
      await this.transport.close();
      console.log('\n✓ Disconnected from DeepWiki MCP server');
    }
  }

  async askDeepWiki(question: string): Promise<string> {
    if (!this.transport)
      throw new Error('Not connected to MCP server. Call connect() first.');


    const q = `You are a helpful assistant with access to documentation for the ${this.repoName} GitHub repository.
      
      VERY IMPORTANT TO REMEMER!!!: you are a UI side assistant for Datagrok platform users, so when asked about how to do something in Datagrok, try to answer with UI in mind and also the code.

      Make sure the output is nicely formatted with markdown syntax where applicable.

      Do not repeat the question in your answer, and do not include the deepwiki links!!!.

      user: ${question}
      `;
    // just call the ask tool directly
    const toolArgs = {
      repoName: this.repoName,
      repository: this.repoName,
      question: q,
    };

    const result = await this.mcpClient.callTool({
      name: this.askTool.name,
      arguments: toolArgs,
    });

    const contentArray = Array.isArray(result.content) ? result.content : [];
    const toolResult = contentArray
      .filter((item: any) => item.type === 'text')
      .map((item: any) => item.text || '')
      .join('\n');
    if (!toolResult || toolResult.length === 0)
      return 'No response generated.';

    // remove deepwiki links from the answer
    const lines = toolResult.split('\n').filter((line) => !line.includes('deepwiki.com')).join('\n');
    return lines;
  }
}

// Example usage
export async function askDeepWiki(question: string) {
  const client = DeepWikiMCPClient.getInstance('datagrok-ai/public');

  try {
    await client.connect();

    const response = await client.askDeepWiki(
      question
    );

    return response;
  } catch (error) {
    console.error('Error:', error);
    return `Error: ${error}`;
  }
}


