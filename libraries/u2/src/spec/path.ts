/* The A4 binding-path grammar (DD5): `$.first`, then `.ident` or `['arbitrary string']` steps —
   dots for identifiers, single-quoted brackets with `\'` and `\\` escapes for everything else.
   No wildcards, no filters, no expressions. Platform-free; the resolver, the editor's validation
   and the binding picker all speak through these three functions. */
import type {Spec, SpecNode} from './spec.js';

const IDENTIFIER = /^[A-Za-z_$][A-Za-z0-9_$]*$/;
const IDENTIFIER_CHAR = /[A-Za-z0-9_$]/;

/** Segments of a well-formed path, the first included; null on anything malformed. */
export function parsePath(path: string): string[] | null {
  if (!path.startsWith('$.'))
    return null;
  const segments: string[] = [];
  let at = 2;
  const ident = (): boolean => {
    let end = at;
    while (end < path.length && IDENTIFIER_CHAR.test(path[end]))
      end++;
    if (end === at || !IDENTIFIER.test(path.slice(at, end)))
      return false;
    segments.push(path.slice(at, end));
    at = end;
    return true;
  };
  if (!ident())
    return null;
  while (at < path.length) {
    if (path[at] === '.') {
      at++;
      if (!ident())
        return null;
    } else if (path[at] === '[') {
      if (path[at + 1] !== '\'')
        return null;
      at += 2;
      let value = '';
      for (;;) {
        if (at >= path.length)
          return null;
        const c = path[at];
        if (c === '\\') {
          const next = path[at + 1];
          if (next !== '\'' && next !== '\\')
            return null;
          value += next;
          at += 2;
        } else if (c === '\'') {
          at++;
          break;
        } else {
          value += c;
          at++;
        }
      }
      if (path[at] !== ']')
        return null;
      at++;
      segments.push(value);
    } else
      return null;
  }
  return segments;
}

/** The picker's assembly: identifier segments dotted, anything else bracketed and escaped —
 * the exact inverse of {@link parsePath}. */
export function segmentsToPath(segments: string[]): string {
  let path = `$.${segments[0]}`;
  for (const segment of segments.slice(1)) {
    path += IDENTIFIER.test(segment) ? `.${segment}` :
      `['${segment.replace(/\\/g, '\\\\').replace(/'/g, '\\\'')}']`;
  }
  return path;
}

/** Whether `path`'s first segment is `name`: `$.old` and `$.old.x` reference it, `$.older` does not. */
export function referencesName(path: string, name: string): boolean {
  const head = `$.${name}`;
  if (!path.startsWith(head))
    return false;
  const rest = path.slice(head.length);
  return rest === '' || rest[0] === '.' || rest[0] === '[';
}

/** Every node — in the form or on the tray — whose binds or events reference `name`: `$.name…`
 * paths, event args included, and `cmd:name.…` component functions. The read-only sibling of
 * `renameRefs`: what captured `name`'s signals at its own render, and must be rendered again
 * when what the name denotes is rebuilt. */
export function referencesOf(spec: Spec, name: string): SpecNode[] {
  const found: SpecNode[] = [];
  const walk = (node: SpecNode): void => {
    if (references(node, name))
      found.push(node);
    for (const child of node.children ?? [])
      walk(child);
  };
  walk(spec.root);
  for (const component of spec.components ?? [])
    walk(component);
  return found;
}

function references(node: SpecNode, name: string): boolean {
  for (const path of Object.values(node.bind ?? {})) {
    if (referencesName(path, name))
      return true;
  }
  for (const entry of Object.values(node.on ?? {})) {
    const cmd = typeof entry === 'string' ? entry : entry?.cmd;
    if (typeof cmd === 'string' && cmd.startsWith(`cmd:${name}.`))
      return true;
    if (typeof entry === 'string')
      continue;
    for (const value of Object.values(entry?.args ?? {})) {
      if (typeof value === 'string' && referencesName(value, name))
        return true;
    }
  }
  return false;
}

/** Bind-time cycle detection (C3): an edge is a named node or component whose bind path starts at
 * another name in the spec — context data cannot cycle. Answers the first loop found, closed
 * (`orders → daysInput → orders`), or null. The editor refuses a patch that would create one;
 * the renderer marks the members of a hand-written one broken. */
export function findBindingCycle(spec: Spec): string[] | null {
  const nodes = new Map<string, SpecNode>();
  const collect = (node: SpecNode): void => {
    if (node.name !== undefined && !nodes.has(node.name))
      nodes.set(node.name, node);
    for (const child of node.children ?? [])
      collect(child);
  };
  collect(spec.root);
  for (const component of spec.components ?? [])
    collect(component);

  const edges = new Map<string, string[]>();
  for (const [name, node] of nodes) {
    const to: string[] = [];
    for (const path of Object.values(node.bind ?? {})) {
      const first = parsePath(path)?.[0];
      if (first !== undefined && nodes.has(first) && !to.includes(first))
        to.push(first);
    }
    edges.set(name, to);
  }

  const state = new Map<string, 'visiting' | 'done'>();
  const stack: string[] = [];
  const visit = (name: string): string[] | null => {
    state.set(name, 'visiting');
    stack.push(name);
    for (const next of edges.get(name) ?? []) {
      if (state.get(next) === 'visiting')
        return [...stack.slice(stack.indexOf(next)), next];
      if (!state.has(next)) {
        const loop = visit(next);
        if (loop)
          return loop;
      }
    }
    stack.pop();
    state.set(name, 'done');
    return null;
  };
  for (const name of nodes.keys()) {
    if (!state.has(name)) {
      const loop = visit(name);
      if (loop)
        return loop;
    }
  }
  return null;
}
