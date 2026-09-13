/// YAML frontmatter: a file whose first line is `---`, up to the next line that is exactly `---`.
import * as yaml from 'js-yaml';

export interface Frontmatter {
  /** Null when the file has no frontmatter or it does not parse (see [error]). */
  data: Record<string, unknown> | null;
  body: string;
  /** 1-based line the body starts at. */
  bodyLine: number;
  /** The raw YAML text, for locating keys by line. */
  yaml: string;
  error?: string;
  /** 1-based line of a YAML error, when known. */
  errorLine?: number;
}

export function splitFrontmatter(text: string): Frontmatter {
  const source = text.charCodeAt(0) === 0xFEFF ? text.slice(1) : text;
  const lines = source.split(/\r?\n/);
  if (lines[0] !== '---')
    return {data: null, body: source, bodyLine: 1, yaml: ''};
  const end = lines.indexOf('---', 1);
  if (end < 0)
    return {data: null, body: source, bodyLine: 1, yaml: '', error: 'frontmatter is not closed by a --- line', errorLine: 1};
  const block = lines.slice(1, end).join('\n');
  const body = lines.slice(end + 1).join('\n');
  try {
    // js-yaml's load rejects duplicate keys; keep it that way (a duplicate key is a silent overwrite otherwise)
    const data = yaml.load(block);
    if (data === null || data === undefined)
      return {data: {}, body, bodyLine: end + 2, yaml: block};
    if (typeof data !== 'object' || Array.isArray(data))
      return {data: null, body, bodyLine: end + 2, yaml: block, error: 'frontmatter must be a YAML mapping', errorLine: 2};
    return {data: data as Record<string, unknown>, body, bodyLine: end + 2, yaml: block};
  } catch (e: any) {
    const line = e.mark ? e.mark.line + 2 : undefined;
    return {data: null, body, bodyLine: end + 2, yaml: block, error: `YAML error: ${e.reason ?? e.message}`, errorLine: line};
  }
}

/** 1-based line of a top-level [key] inside the frontmatter, or undefined. */
export function keyLine(fm: Frontmatter, key: string): number | undefined {
  const lines = fm.yaml.split('\n');
  const index = lines.findIndex((l) => l.startsWith(`${key}:`) || l.startsWith(`"${key}":`) || l.startsWith(`'${key}':`));
  return index < 0 ? undefined : index + 2;
}
