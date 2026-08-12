/** The `.flow` script-entity body: a Datagrok annotation header + the flow JSON document.
 *  Single-writer invariant: only `flowScriptText` produces this text, so header and JSON never disagree. */

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

/** `extras` merges document fields (e.g. `outputViews` layouts) that live outside the graph. */
export function flowScriptText(flow: FlowEditor, settings: FlowSettings,
  extras?: Partial<FuncFlowDocument>): string {
  const tags = settings.tags.includes(FLOW_TAG) ? settings.tags : [...settings.tags, FLOW_TAG];
  const doc = serializeFlow(flow, {...settings, tags});
  if (extras) {
    for (const [k, v] of Object.entries(extras)) {
      if (v !== undefined) (doc as unknown as Record<string, unknown>)[k] = v;
    }
  }
  const header = emitHeaderLines(flow, {
    name: settings.scriptName,
    description: settings.scriptDescription,
    tags,
  }, FLOW_LANGUAGE);
  return header.join('\n') + '\n' + JSON.stringify(doc, null, 2) + '\n';
}

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

export function isFlowBody(text: string | null | undefined): boolean {
  if (!text) return false;
  try {
    parseFlowBody(text);
    return true;
  } catch {
    return false;
  }
}
