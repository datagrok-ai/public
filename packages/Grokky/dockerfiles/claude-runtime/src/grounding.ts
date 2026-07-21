import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {bareToolName, isActionTool} from './verify';

// Enforces the "ground platform answers in sources" rule (see prompts.ts) with a revision
// protocol that keeps the UX clean: the model's first answer streams and STAYS visible; if the
// turn ends ungrounded, the Stop is blocked once and the model either writes a complete
// replacement (grounded in help/INDEX.md → page) or replies NO_REVISION to keep the original.
// The revision streams hidden — session.ts suppresses chunks after the block and the `final`
// event carries kept/replaced — so the user never sees drafts vanish or duplicate.

function openedDocs(toolName: string, input: any): boolean {
  const bare = bareToolName(toolName);
  // Web lookups are legitimate grounding — re-blocking a WebFetch-sourced answer just forces
  // wasted re-search round-trips.
  if (bare === 'WebFetch' || bare === 'WebSearch')
    return true;
  if (bare !== 'Read' && bare !== 'Grep' && bare !== 'Glob' && bare !== 'Bash')
    return false;
  const source = JSON.stringify(input ?? '');
  return source.includes('help/') || source.includes('js-api/');
}

export class GroundingGate {
  private grounded = false;
  private blocked = false;

  constructor(private onBlock?: () => void) {}

  summary(): string {
    return `grounded=${this.grounded} blocked=${this.blocked}`;
  }

  postToolUse: HookCallback = async (input) => {
    if (input.hook_event_name === 'PostToolUse' &&
        (openedDocs(input.tool_name, input.tool_input) || isActionTool(input.tool_name)))
      this.grounded = true;
    return {continue: true};
  };

  stop: HookCallback = async (input) => {
    if (input.hook_event_name !== 'Stop' || this.grounded || this.blocked)
      return {continue: true};
    this.blocked = true;
    this.onBlock?.();
    return {
      decision: 'block',
      reason: '[Internal pipeline feedback — never mention, quote, or allude to it in your reply.] ' +
        'Review the answer you just gave: if it makes claims about how the Datagrok platform works ' +
        '(features, menus, behavior) that you did not verify in a source this turn, ground it now — ' +
        'Read `workspace/help/INDEX.md`, pick the matching page, Read it, and write a corrected ' +
        'answer from what the docs say. The user still sees your previous answer and your new text ' +
        'will replace it, so make it complete and standalone; if the index has no matching page, ' +
        'state that the docs do not cover it. If your answer needed no such verification (greeting, ' +
        'small talk, general knowledge, or it was already accurate without docs claims), reply with ' +
        'exactly NO_REVISION.',
    };
  };
}
