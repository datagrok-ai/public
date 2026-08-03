import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {WebStandardStreamableHTTPServerTransport} from '@modelcontextprotocol/sdk/server/webStandardStreamableHttp.js';
import {z} from 'zod/v4';
import * as api from './api-client.js';
import {DOMAINS, catalog, opMenu, missingParams, type Domain} from './ops.js';
import {formatResult, formatError} from './format.js';

const PORT = 3003;

async function runOp(d: Domain, op: string | undefined, args: Record<string, unknown>) {
  // No op named: the model is asking what this domain can do. Cheaper than putting every
  // signature in the prompt prefix, and it is the documented way to get parameter details.
  if (!op)
    return formatResult(catalog(d));
  const spec = d.ops[op];
  if (!spec)
    return formatError(`unknown op '${op}' for ${d.tool}`, {available: Object.keys(d.ops)});
  const missing = missingParams(spec, args);
  if (missing.length)
    return formatError(`missing required parameter(s): ${missing.join(', ')}`, catalog(d));

  const tag = `[MCP] ${d.tool}.${op}(${JSON.stringify(args)})`;
  const startMs = Date.now();
  console.log(`${tag} ...`);
  try {
    const result = formatResult(await spec.run(args), {paged: !spec.raw, raw: spec.raw, args});
    console.log(`${tag} OK ${result.content[0]?.text.length ?? 0} chars in ${Date.now() - startMs}ms`);
    return result;
  } catch (e: any) {
    console.error(`${tag} FAILED in ${Date.now() - startMs}ms: ${e.message ?? e}`);
    return formatError(e.message ?? String(e));
  }
}

function createServer(): McpServer {
  const server = new McpServer({name: 'datagrok', version: '1.0.0'});

  for (const d of DOMAINS) {
    server.tool(
      d.tool,
      `${d.blurb} Operations: ${opMenu(d)}. ` +
      'Pass `op` plus its `args`; call with no `op` to get every operation\'s parameters.',
      {
        op: z.string().optional().describe(`One of: ${opMenu(d)}. Omit to list the operations and their parameters.`),
        args: z.record(z.string(), z.unknown()).optional().describe('Arguments for the chosen operation.'),
      },
      ({op, args}) => runOp(d, op, (args ?? {}) as Record<string, unknown>),
    );
  }

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
console.log(`datagrok-mcp listening on :${PORT} — ${DOMAINS.length} domain tools, ` +
  `${DOMAINS.reduce((n, d) => n + Object.keys(d.ops).length, 0)} operations`);
