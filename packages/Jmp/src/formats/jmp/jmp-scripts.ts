// Saved table scripts in a .jmp file. Not in JMPReader.jl; the layout below is inferred from the
// "Source" script of the compact_* fixtures and is provisional until more files are read (JMP_PLUGIN.md §4.5).
//
// The scripts sit between the build string and the column data. Each one seen so far is
//   0x03, u32 length, UTF-8 JSL text of that length
// and the text is `Name(body)`, the script's name called on the script.

import {JslNode, JslSyntaxError, isCall, jslKey, parseJsl} from './jsl';
import {JmpInfo} from './jmp-types';

export interface JmpScript {
  /** The script's name in JMP's table panel, e.g. "Source" or "Distribution". */
  name: string;
  /** The script, without the name wrapper. */
  text: string;
  /** The parsed script; undefined when it couldn't be read, with `error` set. */
  body?: JslNode;
  error?: string;
}

const SCRIPT_TAG = 0x03;
const strictUtf8 = new TextDecoder('utf-8', {fatal: true});

/** Finds the table scripts saved in the file. Returns an empty list when there are none or none are recognized. */
export function readTableScripts(bytes: Uint8Array, info: JmpInfo): JmpScript[] {
  const [start, end] = info.propertiesRange;
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const out: JmpScript[] = [];
  for (let p = start; p + 5 < end;) {
    if (bytes[p] === SCRIPT_TAG) {
      const len = view.getUint32(p + 1, true);
      const script = len >= 3 && p + 5 + len <= end ? tryScript(bytes.subarray(p + 5, p + 5 + len)) : null;
      if (script) {
        out.push(script);
        p += 5 + len;
        continue;
      }
    }
    p++;
  }
  return out;
}

function tryScript(raw: Uint8Array): JmpScript | null {
  // cheap checks first: JSL text starts with a name and ends with a parenthesis
  if (!/[A-Za-z_"]/.test(String.fromCharCode(raw[0])) || raw[raw.length - 1] !== 0x29)
    return null;
  let text: string;
  try {
    text = strictUtf8.decode(raw);
  } catch {
    return null;
  }
  if (/[\0-\x08\x0e-\x1f]/.test(text))
    return null;
  // a name JSL can't write bare is saved as `Name("1a fit")(body)`
  const named = /^\s*Name\(\s*"([^"]*)"\s*\)\s*\(/i.exec(text);
  const open = named ? named[0].length - 1 : text.indexOf('(');
  const head = named ? named[1] : text.slice(0, open).trim();
  if (open <= 0 || (!named && !/^("[^"]*"n|[A-Za-z_][\w ]*)$/.test(head)))
    return null;
  const inner = text.slice(open + 1, -1).trim();
  try {
    const root = parseJsl(text);
    if (root.type !== 'call')
      return null;
    if (isCall(root, 'New Script') && root.args[0]?.type === 'string')
      return {name: root.args[0].value, text: textAfterFirstArg(inner), body: root.args[1]};
    const body: JslNode = root.args.length === 1 ? root.args[0] : {type: 'list', items: root.args};
    return {name: root.name, text: inner, body};
  } catch (e) {
    // it has the shape of a script but we can't read it: keep it, so it can be listed
    if (e instanceof JslSyntaxError)
      return {name: head.replace(/^"(.*)"n$/, '$1'), text: inner, error: e.message};
    throw e;
  }
}

function textAfterFirstArg(inner: string): string {
  const m = /^"(?:[^"\\]|\\!.)*"\s*,/.exec(inner);
  return m ? inner.slice(m[0].length).trim() : inner;
}

/** True for a script that only re-imports the table's source file, which JMP saves as "Source". */
export function isSourceScript(s: JmpScript): boolean {
  return jslKey(s.name) === 'source';
}
