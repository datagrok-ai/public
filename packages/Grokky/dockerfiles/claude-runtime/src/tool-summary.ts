import {bareToolName, DOMAIN_TOOLS} from './verify';

interface ToolInputs {
  read_file: {path?: string};
  write_file: {path?: string};
  patch: {path?: string};
  search_files: {pattern?: string; target?: 'content' | 'files'; path?: string};
  terminal: {command?: string};
  process_manage: {action?: string};
  skills_list: {category?: string};
  skill_view: {name?: string};
  web_search: {query?: string};
  web_extract: {urls?: string[]};
  datagrok_exec: {code?: string};
  datagrok_verify: {assertion?: string; description?: string};
  datagrok_show_entities: {entities?: unknown[]};
  AskUserQuestion: {questions?: {question?: string}[]};
}

type ToolName = keyof ToolInputs;

/** The Datagrok MCP server has one tool per domain, each dispatching on `op` (mcp-server/src/ops.ts). */
interface DomainToolInput {
  op?: string;
  args?: Record<string, any>;
}

const inPath = (path?: string) => path ? ` in ${path}` : '';

const toolFormatters: {[K in ToolName]: (i: ToolInputs[K]) => string} = {
  read_file: (i) => `Read ${i.path ?? ''}`,
  write_file: (i) => `Write ${i.path ?? ''}`,
  patch: (i) => `Edit ${i.path ?? ''}`,
  search_files: (i) => i.target === 'files' ? `Find files ${i.pattern ?? ''}${inPath(i.path)}` :
    `Search '${i.pattern ?? ''}'${inPath(i.path)}`,
  terminal: () => 'Run a command',
  process_manage: (i) => `Background process: ${i.action ?? ''}`,
  skills_list: () => 'List skills',
  skill_view: (i) => `Load skill ${i.name ?? ''}`,
  web_search: (i) => `Search the web: ${i.query ?? ''}`,
  web_extract: (i) => `Fetch ${(i.urls ?? []).join(', ')}`,
  datagrok_exec: () => 'Execute in browser',
  datagrok_verify: (i) => `Verify: ${i.description ?? 'action outcome'}`,
  datagrok_show_entities: (i) => {
    const n = (i.entities ?? []).length;
    return `Show ${n} entit${n === 1 ? 'y' : 'ies'}`;
  },
  AskUserQuestion: (i) => `Asking: ${i.questions?.[0]?.question ?? 'a question'}`,
};

function domainSummary(domain: string, i: DomainToolInput): string {
  if (!i.op)
    return `List ${domain} operations`;
  const subject = i.args?.name ?? i.args?.id ?? i.args?.path ?? i.args?.filter ?? '';
  return `${domain}: ${i.op.replace(/_/g, ' ')}${subject ? ' ' + subject : ''}`;
}

/** One line for the panel while a tool runs. View tools and unknown tools read as their name. */
export function toolSummary(name: string, input: Record<string, unknown>): string {
  // Tool names come from the browser and MCP servers: `in` would also find `constructor` or `__proto__`.
  if (Object.prototype.hasOwnProperty.call(toolFormatters, name))
    return toolFormatters[name as ToolName](input as ToolInputs[ToolName]);
  const bare = bareToolName(name);
  return DOMAIN_TOOLS.has(bare) ? domainSummary(bare.replace('datagrok_', ''), input as DomainToolInput) :
    bare.replace(/_/g, ' ');
}
