/* The Grokky AI panel: the platform docks it when `grok.shell.windows.showAI` is on, and Grokky fills
   it with the chat (prompt box, Stop button, responses, AskUserQuestion forms, system messages). */
import {element} from '@datagrok-libraries/bdd';

export const AI_PANEL = '.d4-ai-generation-panel';
export const AI_PROMPT_BOX = `${AI_PANEL} textarea.d4-ai-input-textarea`;

element('AI stop button', {selector: `${AI_PANEL} .d4-ai-panel-run-controls .fa-stop`,
  description: 'the Stop icon the panel shows while a turn runs'});
