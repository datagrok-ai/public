import {OpenAIHelpClient} from '../llm-utils/openAI-client';
import {_package} from '../package';

export interface MCPTool {
  type: 'mcp';
  server_label: string;
  server_description: string;
  server_url: string;
  require_approval?: 'never';
  authorization?: string;
  allowed_tools?: string[];
}

export class MCPServer {
  public readonly name: string;
  public readonly description: string;
  public readonly endpointUrl: string;
  public readonly domain: string;
  public readonly authorization?: string;
  public readonly allowedTools?: string[];
  public readonly requireApproval?: 'never';

  constructor(opts: {
    name: string;
    description: string;
    endpointUrl: string;
    domain: string;
    authorization?: string;
    allowedTools?: string[];
    requireApproval?: 'never';
  }) {
    this.name = opts.name;
    this.description = opts.description;
    this.endpointUrl = opts.endpointUrl;
    this.domain = opts.domain;
    this.authorization = opts.authorization;
    this.allowedTools = opts.allowedTools;
    this.requireApproval = opts.requireApproval;
  }

  async validateConnection(timeoutMs = 5000): Promise<boolean> {
    try {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), timeoutMs);

      const res = await fetch(this.endpointUrl, {method: 'GET', signal: controller.signal});
      clearTimeout(timeout);
      return res.ok;
    } catch (err) {
      _package.logger.error(err);
      return false;
    }
  }
}

export class MCPRegistry {
  private mcpServers: Map<string, MCPServer[]> = new Map();
  private llmClient: OpenAIHelpClient = OpenAIHelpClient.getInstance();

  constructor(opts: { llmClient: OpenAIHelpClient }) {
    if (!opts.llmClient) throw new Error('LLM client is required');
    this.llmClient = opts.llmClient;
  }

  async register(server: MCPServer, timeoutMs = 5000): Promise<void> {
    const isAlive = await server.validateConnection(timeoutMs);
    if (!isAlive) throw new Error(`Unable to connect to MCP server at '${server.endpointUrl}'.`);

    const existing = this.mcpServers.get(server.domain) || [];
    if (existing.some((s) => s.name === server.name))
      throw new Error(`MCP server with name '${server.name}' already registered in domain '${server.domain}'.`);


    this.mcpServers.set(server.domain, [...existing, server]);
  }

  getByDomain(domain: string): MCPServer[] {
    return this.mcpServers.get(domain) || [];
  }

  /**
   * Picks MCP tools for the given prompt based on LLM classification.
   * Returns array of MCPTool objects ready for LLM usage.
   */
  async pickToolsForPrompt(prompt: string, fallbackDomain?: string): Promise<MCPTool[]> {
    let domains: string[] = [];

    try {
      const completion = await this.llmClient.openai.chat.completions.create({
        model: 'gpt-4o-mini',
        messages: [
          {
            role: 'system',
            content: `Classify the user's prompt into ALL applicable domains from the following list:
            \n${[...this.mcpServers.keys()].join(', ')}\n
            Respond only with a JSON array of domain names.`,
          },
          {role: 'user', content: prompt},
        ],
      });

      const responseText = completion.choices[0].message.content?.trim() || '';
      try {
        domains = Array.isArray(JSON.parse(responseText)) ? JSON.parse(responseText) : [responseText];
      } catch {
        domains = [responseText];
      }
    } catch (err) {
      _package.logger.error(err);
      domains = fallbackDomain ? [fallbackDomain] : ['misc'];
    }

    const seen = new Set<string>();
    const tools: MCPTool[] = [];

    for (const domain of domains) {
      for (const server of this.getByDomain(domain)) {
        if (seen.has(server.name)) continue;
        seen.add(server.name);
        const tool: MCPTool = {
          type: 'mcp',
          server_label: server.name,
          server_description: server.description,
          server_url: server.endpointUrl,
        };
        if (server.requireApproval) tool.require_approval = server.requireApproval;
        if (server.authorization) tool.authorization = server.authorization;
        if (server.allowedTools) tool.allowed_tools = server.allowedTools;
        tools.push(tool);
      }
    }

    return tools;
  }
}
