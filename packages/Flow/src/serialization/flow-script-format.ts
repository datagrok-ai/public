/** The `.flow` script-entity body format: a standard Datagrok annotation
 *  header (`//name:` / `//language: flow` / `//input:` …) followed by the
 *  lossless .ffjson document.
 *
 *  The header lets core parse the entity's name and params even when the Flow
 *  package is not installed; the JSON is the single source of truth for the
 *  graph. **Single-writer invariant**: only `flowScriptText` produces this
 *  text, always regenerating header and JSON together from the same graph, so
 *  they can never disagree. */

import {FlowEditor} from '../rete/flow-editor';
import {emitHeaderLines} from '../compiler/script-emitter';
import {FlowSettings, FuncFlowDocument} from './flow-schema';
import {serializeFlow} from './flow-serializer';

export const FLOW_LANGUAGE = 'flow';
export const FLOW_TAG = 'flow';

export interface ParsedFlowBody {
  /** The leading `//`-comment block, verbatim. */
  header: string;
  doc: FuncFlowDocument;
}

/** Serialize the live graph into the `.flow` entity body. */
export function flowScriptText(flow: FlowEditor, settings: FlowSettings): string {
  const tags = settings.tags.includes(FLOW_TAG) ? settings.tags : [...settings.tags, FLOW_TAG];
  const doc = serializeFlow(flow, {...settings, tags});
  const header = emitHeaderLines(flow, {
    name: settings.scriptName,
    description: settings.scriptDescription,
    tags,
  }, FLOW_LANGUAGE);
  return header.join('\n') + '\n' + JSON.stringify(doc, null, 2) + '\n';
}

/** Split a `.flow` entity body back into its header and ffjson document.
 *  Tolerant of blank lines between the two; throws on a missing or
 *  wrong-version JSON payload. */
export function parseFlowBody(text: string): ParsedFlowBody {
  const lines = text.split('\n');
  let i = 0;
  while (i < lines.length && (lines[i].trimStart().startsWith('//') || lines[i].trim() === ''))
    i++;
  const header = lines.slice(0, i).join('\n').trimEnd();
  const json = lines.slice(i).join('\n').trim();
  if (json === '')
    throw new Error('Flow script body has no JSON payload after the header');
  const doc = JSON.parse(json) as FuncFlowDocument;
  if (doc.version !== '2.0')
    throw new Error(`Unsupported flow version "${doc.version}"; expected 2.0`);
  return {header, doc};
}

/** Whether a script body looks like a `.flow` entity body (used by guards). */
export function isFlowBody(text: string | null | undefined): boolean {
  if (!text) return false;
  try {
    parseFlowBody(text);
    return true;
  } catch {
    return false;
  }
}
