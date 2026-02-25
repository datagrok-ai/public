#!/usr/bin/env node
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {WebStandardStreamableHTTPServerTransport} from '@modelcontextprotocol/sdk/server/webStandardStreamableHttp.js';
import {z} from 'zod/v4';
import * as api from './api-client.js';

const PORT = parseInt(process.env['MCP_PORT'] ?? '3000', 10);

type ToolResult = {content: {type: 'text'; text: string}[]};

function formatResult(data: unknown): ToolResult {
  const text = typeof data === 'string' ? data : JSON.stringify(data, null, 2);
  return {content: [{type: 'text', text}]};
}

async function runTool(
  name: string, args: Record<string, unknown>, fn: () => Promise<unknown>,
): Promise<ToolResult> {
  const startMs = Date.now();
  try {
    const result = formatResult(await fn());
    console.log(`[MCP] ${name}(${JSON.stringify(args)}) OK in ${Date.now() - startMs}ms`);
    return result;
  } catch (e) {
    console.error(`[MCP] ${name}(${JSON.stringify(args)}) FAILED in ${Date.now() - startMs}ms:`, e);
    throw e;
  }
}

function createServer(): McpServer {
  const server = new McpServer({name: 'datagrok', version: '1.0.0'});

  server.tool(
    'list_functions', 'List Datagrok functions with optional smart filter',
    {filter: z.string().optional().describe('Smart search filter (e.g. source="script")')},
    ({filter}) => runTool('list_functions', {filter}, () => api.listFunctions(filter)),
  );

  server.tool(
    'get_function', 'Get detailed information about a function',
    {id: z.string().describe('Function ID or Grok name (e.g. "Namespace:FuncName")')},
    ({id}) => runTool('get_function', {id}, () => api.getFunction(id)),
  );

  server.tool(
    'call_function', 'Execute a Datagrok function',
    {
      name: z.string().describe('Grok name of the function (e.g. "Namespace:FuncName")'),
      params: z.record(z.string(), z.unknown()).optional()
        .describe('Parameters to pass to the function'),
    },
    ({name, params}) => runTool('call_function', {name, params},
      () => api.callFunction(name, params)),
  );

  server.tool(
    'list_scripts', 'List all scripts', {},
    () => runTool('list_scripts', {}, () => api.listFunctions('source="script"')),
  );

  server.tool(
    'list_queries', 'List all data queries', {},
    () => runTool('list_queries', {}, () => api.listFunctions('source="data-query"')),
  );

  server.tool(
    'create_script', 'Create a new script in Datagrok',
    {
      script: z.string().describe('Script source code'),
      name: z.string().optional().describe('Script name'),
      language: z.string().optional()
        .describe('Script language (python, r, julia, nodejs, grok, octave). Defaults to python'),
    },
    ({script, name, language}) => runTool('create_script', {name, language}, () => api.saveFunction({
      source: 'script', script,
      name: name ?? 'Untitled Script',
      language: language ?? 'python',
    })),
  );

  server.tool(
    'create_query', 'Create a new data query in Datagrok',
    {
      connectionId: z.string().describe('ID of the data connection'),
      query: z.string().describe('SQL query text'),
      name: z.string().optional().describe('Query name'),
    },
    ({connectionId, query, name}) => runTool('create_query', {connectionId, name}, () => api.saveFunction({
      source: 'data-query', query,
      connection: {id: connectionId},
      name: name ?? 'Untitled Query',
    })),
  );

  server.tool(
    'list_files', 'List files in a connector directory',
    {
      connector: z.string().describe('Connector Grok name (e.g. "System:DemoFiles")'),
      path: z.string().optional().describe('Directory path within the connector'),
    },
    ({connector, path}) => runTool('list_files', {connector, path},
      () => api.listFiles(connector, path)),
  );

  server.tool(
    'download_file', 'Download a file from a connector',
    {
      connector: z.string().describe('Connector Grok name (e.g. "System:DemoFiles")'),
      path: z.string().describe('File path within the connector'),
    },
    ({connector, path}) => runTool('download_file', {connector, path},
      () => api.downloadFile(connector, path)),
  );

  server.tool(
    'upload_file', 'Upload content to a file in a connector',
    {
      connector: z.string().describe('Connector Grok name (e.g. "System:DemoFiles")'),
      path: z.string().describe('File path within the connector'),
      content: z.string().describe('File content to upload'),
    },
    ({connector, path, content}) => runTool('upload_file', {connector, path},
      () => api.uploadFile(connector, path, content)),
  );

  server.tool(
    'whoami', 'Get current authenticated user info', {},
    () => runTool('whoami', {}, () => api.getCurrentUser()),
  );

  return server;
}

const app = new Hono();

app.get('/health', (c) => c.json({status: 'ok'}));

app.all('/mcp', async (c) => {
  console.log(`[MCP] ${c.req.method} /mcp`);
  const transport = new WebStandardStreamableHTTPServerTransport();
  const server = createServer();
  await server.connect(transport);
  return transport.handleRequest(c.req.raw);
});

serve({fetch: app.fetch, port: PORT});
console.log(`datagrok-mcp listening on :${PORT}`);
