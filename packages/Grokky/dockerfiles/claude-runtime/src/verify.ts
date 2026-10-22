import type {TurnGate} from './types';

const MAX_VERIFY_BLOCKS = 3;

const READONLY_NAME_RE = /^(whoami$|list|get|search|read|download|find)/;
const READONLY_EXTRAS = new Set(['datagrok_show_entities']);

// The datagrok MCP server dispatches on an `op` argument (one tool per domain), so the tool name
// alone no longer says whether a call mutates anything — `datagrok_spaces` is both `list` and
// `delete`. Read the op when there is one and fall back to the name otherwise. Fail-closed: an
// unrecognised call counts as an action, so the verify gate over-asks rather than letting an
// unverified mutation through.
export const DOMAIN_TOOLS = new Set([
  'datagrok_functions', 'datagrok_files', 'datagrok_projects', 'datagrok_spaces', 'datagrok_platform',
]);

export function isReadonlyTool(bare: string, input?: unknown): boolean {
  if (READONLY_EXTRAS.has(bare))
    return true;
  if (DOMAIN_TOOLS.has(bare)) {
    const op = (input as {op?: unknown} | undefined)?.op;
    // No op is a catalog request — pure discovery, nothing to verify.
    return typeof op === 'string' ? READONLY_NAME_RE.test(op) : true;
  }
  return READONLY_NAME_RE.test(bare);
}

export function bareToolName(name: string): string {
  return name.replace(/^mcp__.+?__/, '');
}

// The agent's own tools change files, not the platform, and AskUserQuestion only asks. Every other tool (browser,
// view and Datagrok MCP tools) can, and counts as an action unless it reads.
const AGENT_TOOLS = new Set(['read_file', 'write_file', 'patch', 'search_files', 'terminal', 'process_manage',
  'skills_list', 'skill_view', 'skill_manage', 'web_search', 'web_extract', 'AskUserQuestion']);

export function isActionTool(toolName: string, input?: unknown): boolean {
  const bare = bareToolName(toolName);
  return !AGENT_TOOLS.has(bare) && !isReadonlyTool(bare, input);
}

function parseResult(text: string): {passed?: unknown; verified?: {passed?: unknown}} | undefined {
  try {
    return JSON.parse(text);
  } catch {
    return undefined;
  }
}

export interface VerifyStats {
  actions: number;
  verifies: number;
  passes: number;
  blocks: number;
  exhausted: boolean;
}

export class Verifier implements TurnGate {
  private pendingActions = 0;
  private blockCount = 0;
  private verifyFailures = 0;
  private stats: VerifyStats = {actions: 0, verifies: 0, passes: 0, blocks: 0, exhausted: false};

  constructor(private onBlock?: () => void) {}

  get exhausted(): boolean {
    return this.stats.exhausted;
  }

  get hadActions(): boolean {
    return this.stats.actions > 0;
  }

  statsLine(): string {
    return JSON.stringify(this.stats);
  }

  afterTool(toolName: string, args: Record<string, unknown>, result: string): void {
    const bare = bareToolName(toolName);
    if (bare === 'datagrok_verify') {
      this.stats.verifies++;
      const res = parseResult(result);
      if (res?.passed === true) {
        this.stats.passes++;
        this.pendingActions = 0;
        this.blockCount = 0;
        this.verifyFailures = 0;
      } else
        this.verifyFailures++;
    } else if (isActionTool(toolName, args)) {
      this.stats.actions++;
      // datagrok_exec can self-verify: the browser runs the provided verify.assertion right after
      // the action code (same round-trip) and returns it as `verified` — a passing one is exactly
      // a datagrok_verify pass, minus the extra model round-trip, so it clears the gate the same way.
      const res = bare === 'datagrok_exec' ? parseResult(result) : undefined;
      if (res?.verified?.passed === true) {
        this.stats.verifies++;
        this.stats.passes++;
        this.pendingActions = 0;
        this.blockCount = 0;
        this.verifyFailures = 0;
      } else
        this.pendingActions++;
    }
  }

  atStop(): string | undefined {
    if (this.pendingActions === 0)
      return undefined;
    if (this.blockCount >= MAX_VERIFY_BLOCKS) {
      this.stats.exhausted = true;
      console.warn(`verify: gate exhausted after ${MAX_VERIFY_BLOCKS} blocks — ` +
        `${this.pendingActions} action(s) end unverified`);
      return undefined;
    }
    this.blockCount++;
    this.stats.blocks++;
    this.onBlock?.();
    const internal = '[Internal pipeline feedback — never mention, quote, or allude to it in your reply.] ' +
      'The user still sees your previous answer. Once a verify passes: if that answer was fully ' +
      'accurate, reply with exactly NO_REVISION; otherwise write a complete standalone replacement ' +
      'reporting the observed result — it replaces the previous answer. ';
    return this.verifyFailures > 0 ?
      internal +
      `Verification has failed ${this.verifyFailures} time(s). A failed verify means EITHER the action ` +
      'did not take effect OR the assertion itself was faulty (most commonly a missing `return`). ' +
      'FIRST re-read the affected state with a read-only check to see what actually happened. ' +
      'NEVER redo a state-changing action unless the re-read shows it is genuinely missing — actions ' +
      'are not idempotent, and redoing one that already worked creates duplicates. Then call ' +
      'datagrok_verify again with a corrected assertion. Do not report success until a verify passes.' :
      internal +
      'You performed one or more actions but have not verified they took effect. Call datagrok_verify ' +
      'with an assertion that RE-READS all affected state from t/view/grok. Do NOT redo the actions ' +
      'themselves. If verification fails, fix the problem — do not report success.';
  }
}
