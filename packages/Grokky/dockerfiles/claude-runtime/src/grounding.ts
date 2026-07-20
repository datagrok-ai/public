import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {bareToolName, isActionTool} from './verify';

function openedDocs(toolName: string, input: any): boolean {
  const bare = bareToolName(toolName);
  // Web lookups are legitimate grounding — re-blocking a WebFetch-sourced answer just forces
  // wasted re-search round-trips (measured: each block costs a full discarded answer).
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
    return {
      decision: 'block',
      reason: 'You answered a platform question without opening the documentation this turn. A skill ' +
        'alone is not enough — Read `workspace/help/INDEX.md`, pick the matching page from it, Read that ' +
        'page, and quote the menu path / syntax from it. If the index has no page for the topic, say the ' +
        'docs do not cover it — do not keep searching. If the user\'s message is not a platform how-to ' +
        'question, disregard this and answer normally.',
    };
  };
}
