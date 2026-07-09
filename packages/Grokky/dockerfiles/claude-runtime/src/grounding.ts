import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {bareToolName, isActionTool} from './verify';

function openedDocs(toolName: string, input: any): boolean {
  const bare = bareToolName(toolName);
  if (bare !== 'Read' && bare !== 'Grep' && bare !== 'Glob' && bare !== 'Bash')
    return false;
  return JSON.stringify(input ?? '').includes('help/');
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
        'alone is not enough — grep `workspace/help/` for the terms in the question, read the page, and ' +
        'quote the menu path / syntax from it. If a help search finds nothing, say the docs do not cover ' +
        'it. If the user\'s message is not a platform how-to question, disregard this and answer normally.',
    };
  };
}
