import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {WebStandardStreamableHTTPServerTransport} from '@modelcontextprotocol/sdk/server/webStandardStreamableHttp.js';
import {z} from 'zod/v4';
import * as api from './api-client.js';

const PORT = 3003;

type ToolResult = {content: {type: 'text'; text: string}[]};

function formatResult(data: unknown): ToolResult {
  const text = typeof data === 'string' ? data : JSON.stringify(data, null, 2);
  return {content: [{type: 'text', text}]};
}

async function runTool(
  name: string, args: Record<string, unknown>, fn: () => Promise<unknown>,
): Promise<ToolResult> {
  const tag = `[MCP] ${name}(${JSON.stringify(args)})`;
  const startMs = Date.now();
  console.log(`${tag} ...`);
  try {
    const result = formatResult(await fn());
    const len = result.content[0]?.text.length ?? 0;
    console.log(`${tag} OK ${len} chars in ${Date.now() - startMs}ms`);
    return result;
  } catch (e: any) {
    console.error(`${tag} FAILED in ${Date.now() - startMs}ms: ${e.message ?? e}`);
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
    {id: z.string().describe('Function ID or name (e.g. "Namespace:FuncName")')},
    ({id}) => runTool('get_function', {id}, () => api.getFunction(id)),
  );

  server.tool(
    'call_function', 'Execute a Datagrok function',
    {
      name: z.string().describe('Name of the function (e.g. "Namespace:FuncName")'),
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
      connector: z.string().describe('Connector name (e.g. "System:DemoFiles")'),
      path: z.string().optional().describe('Directory path within the connector'),
    },
    ({connector, path}) => runTool('list_files', {connector, path},
      () => api.listFiles(connector, path)),
  );

  server.tool(
    'download_file', 'Download a file from a connector',
    {
      connector: z.string().describe('Connector name (e.g. "System:DemoFiles")'),
      path: z.string().describe('File path within the connector'),
    },
    ({connector, path}) => runTool('download_file', {connector, path},
      () => api.downloadFile(connector, path)),
  );

  server.tool(
    'upload_file', 'Upload content to a file in a connector',
    {
      connector: z.string().describe('Connector name (e.g. "System:DemoFiles")'),
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

  // --- Projects ---

  server.tool(
    'list_projects', 'List projects with optional smart filter',
    {filter: z.string().optional().describe('Smart search filter')},
    ({filter}) => runTool('list_projects', {filter}, () => api.listProjects(filter)),
  );

  server.tool(
    'get_project', 'Get detailed information about a project',
    {id: z.string().describe('Project ID')},
    ({id}) => runTool('get_project', {id}, () => api.getProject(id)),
  );

  server.tool(
    'create_project', 'Create a new project',
    {
      name: z.string().describe('Project name'),
      description: z.string().optional().describe('Project description'),
    },
    ({name, description}) => runTool('create_project', {name},
      () => api.createProject(name, description)),
  );

  server.tool(
    'delete_project', 'Delete a project',
    {id: z.string().describe('Project ID')},
    ({id}) => runTool('delete_project', {id}, () => api.deleteProject(id)),
  );

  server.tool(
    'search_project', 'Search for a project by name',
    {name: z.string().describe('Project name to search for')},
    ({name}) => runTool('search_project', {name}, () => api.searchProject(name)),
  );

  server.tool(
    'list_recent_projects', 'List recently accessed projects', {},
    () => runTool('list_recent_projects', {}, () => api.listRecentProjects()),
  );

  // --- Spaces ---

  server.tool(
    'list_spaces', 'List root spaces with optional smart filter',
    {filter: z.string().optional().describe('Smart search filter')},
    ({filter}) => runTool('list_spaces', {filter}, () => api.listSpaces(filter)),
  );

  server.tool(
    'get_space', 'Get detailed information about a space',
    {id: z.string().describe('Space ID')},
    ({id}) => runTool('get_space', {id}, () => api.getSpace(id)),
  );

  server.tool(
    'create_space', 'Create a new root space',
    {name: z.string().describe('Space name')},
    ({name}) => runTool('create_space', {name}, () => api.createRootSpace(name)),
  );

  server.tool(
    'delete_space', 'Delete a space',
    {id: z.string().describe('Space ID')},
    ({id}) => runTool('delete_space', {id}, () => api.deleteSpace(id)),
  );

  server.tool(
    'create_subspace', 'Create a subspace within a space',
    {
      spaceId: z.string().describe('Parent space ID'),
      name: z.string().describe('Subspace name'),
      link: z.boolean().optional().describe('Create as a link reference instead of owned child'),
    },
    ({spaceId, name, link}) => runTool('create_subspace', {spaceId, name},
      () => api.createSubspace(spaceId, name, link)),
  );

  server.tool(
    'list_space_children', 'List children of a space (subspaces, entities, etc.)',
    {
      spaceId: z.string().describe('Space ID'),
      types: z.string().optional()
        .describe('Comma-separated entity types to filter (e.g. "Script,DataQuery,Project")'),
      includeLinked: z.boolean().optional().describe('Include linked (non-owned) children'),
    },
    ({spaceId, types, includeLinked}) => runTool('list_space_children', {spaceId, types},
      () => api.listSpaceChildren(spaceId, types, includeLinked)),
  );

  server.tool(
    'add_entity_to_space', 'Add an entity (script, query, connection, etc.) to a space',
    {
      spaceId: z.string().describe('Space ID'),
      entityId: z.string().describe('Entity ID to add'),
      link: z.boolean().optional().describe('Add as a link reference instead of owned child'),
    },
    ({spaceId, entityId, link}) => runTool('add_entity_to_space', {spaceId, entityId},
      () => api.addEntityToSpace(spaceId, entityId, link)),
  );

  server.tool(
    'remove_entity_from_space', 'Remove an entity from a space',
    {
      spaceId: z.string().describe('Space ID'),
      entityId: z.string().describe('Entity ID to remove'),
    },
    ({spaceId, entityId}) => runTool('remove_entity_from_space', {spaceId, entityId},
      () => api.removeEntityFromSpace(spaceId, entityId)),
  );

  server.tool(
    'read_space_file', 'Read a file from space storage',
    {
      spaceId: z.string().describe('Space ID'),
      path: z.string().describe('File path within the space'),
    },
    ({spaceId, path}) => runTool('read_space_file', {spaceId, path},
      () => api.readSpaceFile(spaceId, path)),
  );

  server.tool(
    'write_space_file', 'Write a file to space storage',
    {
      spaceId: z.string().describe('Space ID'),
      path: z.string().describe('File path within the space'),
      content: z.string().describe('File content to write'),
    },
    ({spaceId, path, content}) => runTool('write_space_file', {spaceId, path},
      () => api.writeSpaceFile(spaceId, path, content)),
  );

  server.tool(
    'delete_space_file', 'Delete a file from space storage',
    {
      spaceId: z.string().describe('Space ID'),
      path: z.string().describe('File path within the space'),
    },
    ({spaceId, path}) => runTool('delete_space_file', {spaceId, path},
      () => api.deleteSpaceFile(spaceId, path)),
  );

  return server;
}

const app = new Hono();

app.get('/health', (c) => c.json({status: 'ok'}));

app.all('/mcp', async (c) => {
  const apiKey = c.req.header('x-user-api-key') ?? '';
  const apiUrl = c.req.header('x-datagrok-api-url') ?? '';
  const handler = async () => {
    const transport = new WebStandardStreamableHTTPServerTransport();
    const server = createServer();
    await server.connect(transport);
    return transport.handleRequest(c.req.raw);
  };
  return api.runWithContext({apiKey, apiUrl}, handler);
});

serve({fetch: app.fetch, port: PORT});
console.log(`datagrok-mcp listening on :${PORT}`);
