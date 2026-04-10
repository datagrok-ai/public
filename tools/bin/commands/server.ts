import {NodeDapi} from '../utils/node-dapi';
import {createClient} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';

const ENTITIES = ['users', 'groups', 'functions', 'connections', 'queries', 'scripts', 'packages', 'reports', 'files'];
const VERBS = ['list', 'get', 'delete'];

export async function server(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1);
  const entity: string | undefined = args[0];
  const verb: string | undefined = args[1];
  const rest: string[] = args.slice(2);

  const output: OutputFormat = argv.output ?? argv.o ?? 'table';
  const limit: number = Number(argv.limit ?? argv.l ?? 50);
  const offset: number = Number(argv.offset ?? 0);
  const filter: string = argv.filter ?? argv.f ?? '';
  const host: string | undefined = argv.host;
  const recursive: boolean = !!(argv.r ?? argv.recursive);

  if (!entity || argv.help) {
    console.log(HELP_SERVER);
    return true;
  }

  let client;
  try {
    client = await createClient(host);
  } catch (err: any) {
    printError(err);
    return false;
  }
  const dapi = new NodeDapi(client);

  try {
    if (entity === 'raw') return handleRaw(dapi, verb, rest, output);
    if (entity === 'describe') return handleDescribe(dapi, verb ?? rest[0], output);
    if (entity === 'functions' && verb === 'run') return handleFuncRun(dapi, rest, argv, output);
    if (entity === 'files' && verb === 'list') {
      const path = rest[0] ?? '';
      const result = await dapi.files.list(path, recursive);
      printOutput(result, output);
      return true;
    }
    if (entity === 'files' && verb === 'get') {
      const result = await dapi.files.get(rest[0]);
      printOutput(result, output);
      return true;
    }
    if (entity === 'files' && verb === 'delete') {
      await dapi.files.delete(rest[0]);
      if (output !== 'quiet') console.log(`Deleted ${rest[0]}`);
      return true;
    }

    const source = (dapi as any)[entity];
    if (!source || !ENTITIES.includes(entity)) {
      printError(new Error(`Unknown entity type: '${entity}'. Valid: ${ENTITIES.join(', ')}`));
      return false;
    }
    if (!verb) {
      console.log(`Usage: grok s ${entity} <${VERBS.join('|')}> [args]`);
      return true;
    }

    if (verb === 'list') {
      const page = Math.floor(offset / limit);
      const results = await source
        .filter(filter)
        .by(limit)
        .page(page)
        .list();
      printOutput(results, output);
      return true;
    }
    if (verb === 'get') {
      if (!rest[0]) { printError(new Error('Usage: grok s <entity> get <id>')); return false; }
      const result = await source.find(rest[0]);
      printOutput(result, output);
      return true;
    }
    if (verb === 'delete') {
      if (!rest[0]) { printError(new Error('Usage: grok s <entity> delete <id>')); return false; }
      await source.delete(rest[0]);
      if (output !== 'quiet') console.log(`Deleted ${entity}/${rest[0]}`);
      return true;
    }

    printError(new Error(`Unknown verb: '${verb}'. Valid: ${VERBS.join(', ')}`));
    return false;
  } catch (err: any) {
    printError(err);
    return false;
  }
}

async function handleRaw(dapi: NodeDapi, method: string | undefined, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!method || !rest[0]) {
    printError(new Error('Usage: grok s raw <METHOD> <path>'));
    return false;
  }
  const path = rest[0];
  const result = await dapi.raw(method, path);
  printOutput(result, output);
  return true;
}

async function handleDescribe(dapi: NodeDapi, entityType: string | undefined, output: OutputFormat): Promise<boolean> {
  if (!entityType) {
    printError(new Error('Usage: grok s describe <entity-type>'));
    return false;
  }
  const result = await dapi.describe(entityType);
  printOutput(result, output);
  return true;
}

async function handleFuncRun(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  let funcName = rest[0];
  if (!funcName) {
    printError(new Error('Usage: grok s functions run <Name:function(args...)> [--json params.json]'));
    return false;
  }

  let params: Record<string, any> = {};

  if (argv.json) {
    const fs = require('fs');
    try {
      params = JSON.parse(fs.readFileSync(argv.json, 'utf8'));
    } catch (err: any) {
      printError(new Error(`Cannot read params file: ${err.message}`));
      return false;
    }
  } else {
    const parsed = parseFuncCall(funcName);
    funcName = parsed.name;
    params = parsed.params;
  }

  const result = await dapi.functions.run(funcName, params);
  printOutput(result, output);
  return true;
}

export function parseFuncCall(expr: string): {name: string; params: Record<string, any>} {
  const parenIdx = expr.indexOf('(');
  if (parenIdx === -1) return {name: expr, params: {}};

  const name = expr.slice(0, parenIdx);
  const inner = expr.slice(parenIdx + 1, expr.lastIndexOf(')')).trim();

  if (!inner) return {name, params: {}};

  if (inner.startsWith('{')) {
    try {
      const params = JSON.parse(inner.replace(/([{,]\s*)([a-zA-Z_]\w*)(\s*:)/g, '$1"$2"$3'));
      return {name, params};
    } catch {
      return {name, params: {}};
    }
  }

  const positional = inner.split(',').map((s) => {
    s = s.trim();
    if ((s.startsWith('"') && s.endsWith('"')) || (s.startsWith("'") && s.endsWith("'")))
      return s.slice(1, -1);
    const n = Number(s);
    return isNaN(n) ? s : n;
  });
  return {name, params: Object.fromEntries(positional.map((v, i) => [String(i), v]))};
}

const HELP_SERVER = `
Usage: grok server <entity> <verb> [args] [options]
       grok s <entity> <verb> [args] [options]

Manage a Datagrok server from the command line.

Entities:
  users, groups, functions, connections, queries, scripts, packages, reports, files

Verbs:
  list      List entities
  get       Get a single entity by ID or name
  delete    Delete an entity by ID

Special commands:
  grok s functions run <Name:func(args)>   Call a function
  grok s files list [path] [-r]            List files (recursive with -r)
  grok s raw <METHOD> <path>               Hit any API endpoint
  grok s describe <entity-type>            Show entity JSON schema

Options:
  --host <alias|url>    Server alias from config or full URL
  --output <format>     Output format: table (default), json, csv, quiet
  --filter <text>       Smart filter expression
  --limit <n>           Page size (default: 50)
  --offset <n>          Start offset (default: 0)
  -r, --recursive       Recursive (for files list)
  --json <file>         Read function parameters from JSON file

Examples:
  grok s users list
  grok s users list --output json --limit 10
  grok s connections list --filter "PostgreSQL"
  grok s connections get <id>
  grok s connections delete <id>
  grok s functions run 'Chem:smilesToMw("ccc")'
  grok s functions run Chem:test --json params.json
  grok s files list "System:AppData" -r
  grok s raw GET /api/users/current
  grok s describe connections
  grok s users list --host dev
  grok s users list --host "https://mygrok.com/api"
`;
