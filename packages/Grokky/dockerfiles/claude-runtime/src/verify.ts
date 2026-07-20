import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';

const MAX_VERIFY_BLOCKS = 3;

const READONLY_NAME_RE = /^(whoami$|list_|get_|search_|read_|download_)/;
const READONLY_EXTRAS = new Set(['datagrok_show_entities']);

export function isReadonlyTool(bare: string): boolean {
  return READONLY_NAME_RE.test(bare) || READONLY_EXTRAS.has(bare);
}

export function bareToolName(name: string): string {
  return name.replace(/^mcp__.+?__/, '');
}

export function isActionTool(toolName: string): boolean {
  const bare = bareToolName(toolName);
  return bare === 'datagrok_exec' || (toolName.startsWith('mcp__') && !isReadonlyTool(bare));
}

function parseMcpToolResponse(resp: unknown): any {
  try {
    const content = (resp as any)?.content ?? resp;
    const text = Array.isArray(content) ? content[0]?.text : undefined;
    return typeof text === 'string' ? JSON.parse(text) : resp;
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

export class Verifier {
  private pendingActions = 0;
  private blockCount = 0;
  private verifyFailures = 0;
  private stats: VerifyStats = {actions: 0, verifies: 0, passes: 0, blocks: 0, exhausted: false};

  get exhausted(): boolean {
    return this.stats.exhausted;
  }

  get hadActions(): boolean {
    return this.stats.actions > 0;
  }

  statsLine(): string {
    return JSON.stringify(this.stats);
  }

  postToolUse: HookCallback = async (input) => {
    if (input.hook_event_name !== 'PostToolUse')
      return {continue: true};
    const bare = bareToolName(input.tool_name);
    if (bare === 'datagrok_verify') {
      this.stats.verifies++;
      const res = parseMcpToolResponse(input.tool_response);
      if (res?.passed === true) {
        this.stats.passes++;
        this.pendingActions = 0;
        this.blockCount = 0;
        this.verifyFailures = 0;
      } else
        this.verifyFailures++;
    } else if (isActionTool(input.tool_name)) {
      this.stats.actions++;
      // datagrok_exec can self-verify: the browser runs the provided verify.assertion right after
      // the action code (same round-trip) and returns it as `verified` — a passing one is exactly
      // a datagrok_verify pass, minus the extra model round-trip.
      const res = bare === 'datagrok_exec' ? parseMcpToolResponse(input.tool_response) : undefined;
      if (res?.verified?.passed === true) {
        this.stats.verifies++;
        this.stats.passes++;
      } else
        this.pendingActions++;
    }
    return {continue: true};
  };

  stop: HookCallback = async (input) => {
    if (input.hook_event_name !== 'Stop')
      return {continue: true};
    if (this.pendingActions === 0)
      return {continue: true};
    if (this.blockCount >= MAX_VERIFY_BLOCKS) {
      this.stats.exhausted = true;
      console.warn(`verify: gate exhausted after ${MAX_VERIFY_BLOCKS} blocks — ` +
        `${this.pendingActions} action(s) end unverified`);
      return {continue: true};
    }
    this.blockCount++;
    this.stats.blocks++;
    return {
      decision: 'block',
      reason: this.verifyFailures > 0 ?
        `Verification has failed ${this.verifyFailures} time(s) — the action did NOT take effect as claimed. ` +
        'Fix the problem, then call datagrok_verify again with an assertion that re-reads the affected ' +
        'state. Do not report success until a verify passes.' :
        'You performed one or more actions but have not verified they took effect. Call datagrok_verify ' +
        'with an assertion that RE-READS all affected state from t/view/grok, then report the observed ' +
        'result. If verification fails, fix the problem — do not report success.',
    };
  };
}
