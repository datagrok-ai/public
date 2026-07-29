import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {bareToolName, isActionTool} from './verify';

// Enforces the "ground platform answers in sources" rule (see prompts.ts) with a revision
// protocol that keeps the UX clean: the model's first answer streams and STAYS visible; if the
// turn ends ungrounded, the Stop is blocked once and the model either writes a complete
// replacement (grounded in help/INDEX.md → page) or replies NO_REVISION to keep the original.
// The revision streams hidden — session.ts suppresses chunks after the block and the `final`
// event carries kept/replaced — so the user never sees drafts vanish or duplicate.

// Pure small talk never needs grounding, and the gate is not free: a Stop-block costs a full
// hidden revision call (~28k-token prefix re-read) that ends in NO_REVISION, and past 1.2s the
// panel shows "Revising…" over a greeting. So the gate is not armed at all for messages that are
// only a greeting/thanks/acknowledgement. The check runs on the USER'S OWN text: the panel
// prepends a workspace-context block terminated by '\n---\n\n' (see panel.ts
// prependViewContext), so with a dashboard open the message is '<context>\n---\n\nhello' —
// strip the context first or nothing ever matches. Conservative on purpose: anything beyond
// small talk ("hi, how do I add a scatter plot?") must keep the gate.
const SMALL_TALK_RE = new RegExp(
  '^(?:hi|hii+|hello|hey|yo|good\\s+(?:morning|afternoon|evening)|how\\s+are\\s+you|' +
  'thanks?|thank\\s+you|thx|ty|ok(?:ay)?|cool|nice|great|awesome|perfect|got\\s+it|' +
  'bye|goodbye|see\\s+you|good\\s+night)' +
  '(?:\\s+(?:there|again|grokky|claude))?[\\s.!?,]*$', 'i');

// The context block always opens with this exact header (exec-blocks.ts buildWorkspaceContext),
// so stripping is anchored on it — a user message that merely contains '---' is never split,
// which keeps misclassification fail-closed (gate stays armed).
const CONTEXT_HEADER = 'Workspace state (live, changes as the user navigates):';

export function isSmallTalk(message: string): boolean {
  let text = message;
  if (text.startsWith(CONTEXT_HEADER)) {
    const sep = text.indexOf('\n---\n\n');
    if (sep < 0)
      return false;
    text = text.slice(sep + 6);
  }
  text = text.trim();
  return text.length <= 60 && SMALL_TALK_RE.test(text);
}

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
        (openedDocs(input.tool_name, input.tool_input) || isActionTool(input.tool_name, input.tool_input)))
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
