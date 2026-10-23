// A reader for the declarative part of JSL (JMP Scripting Language), enough to read a saved table script
// such as `Graph Builder(Variables(X(:Age), Y(:Weight)), Elements(Points(X, Y)))` into a tree.
// It never evaluates anything. See JMP_PLUGIN.md §11 (Phase 3).

export class JslSyntaxError extends Error {
  constructor(message: string, readonly pos: number) {
    super(`${message} (at character ${pos + 1})`);
    this.name = 'JslSyntaxError';
  }
}

export type JslNode =
  | {type: 'call', name: string, args: JslNode[]}
  | {type: 'name', name: string}
  | {type: 'column', name: string}
  | {type: 'string', value: string}
  | {type: 'number', value: number}
  | {type: 'list', items: JslNode[]}
  | {type: 'index', target: JslNode, args: JslNode[]}
  | {type: 'unary', op: string, operand: JslNode}
  | {type: 'binary', op: string, left: JslNode, right: JslNode};

/** JSL ignores case and whitespace in names: `Graph Builder`, `graphbuilder` and `GRAPH BUILDER` are one name. */
export function jslKey(name: string): string {
  return name.replace(/\s+/g, '').toLowerCase();
}

type Token =
  | {kind: 'name', text: string, pos: number}
  | {kind: 'string', text: string, pos: number}
  | {kind: 'nameLiteral', text: string, pos: number}
  | {kind: 'number', text: string, pos: number}
  | {kind: 'punct', text: string, pos: number}
  | {kind: 'end', text: '', pos: number};

const PUNCT = ['<<', '>>', '==', '!=', '<=', '>=', '||', '&&', '::', '+=', '-=', '*=', '/=',
  '(', ')', '{', '}', '[', ']', ',', ';', ':', '=', '<', '>', '+', '-', '*', '/', '!', '|', '&', '^', '?'];

const ESCAPES: Record<string, string> = {'"': '"', 'n': '\n', 'r': '\r', 't': '\t', '\\': '\\', '0': '\0', 'N': '\n'};

function tokenize(src: string): Token[] {
  const tokens: Token[] = [];
  let i = 0;
  while (i < src.length) {
    const c = src[i];
    if (/\s/.test(c)) {
      i++;
      continue;
    }
    if (src.startsWith('//', i)) {
      while (i < src.length && src[i] !== '\n')
        i++;
      continue;
    }
    if (src.startsWith('/*', i)) {
      const end = src.indexOf('*/', i + 2);
      if (end < 0)
        throw new JslSyntaxError('Unclosed comment', i);
      i = end + 2;
      continue;
    }
    const pos = i;
    if (src.startsWith('"\\[', i)) {
      // raw string: "\[ ... ]\"
      const end = src.indexOf(']\\"', i + 3);
      if (end < 0)
        throw new JslSyntaxError('Unclosed string', pos);
      tokens.push({kind: 'string', text: src.slice(i + 3, end), pos});
      i = end + 3;
      continue;
    }
    if (c === '"') {
      let text = '';
      i++;
      for (;;) {
        if (i >= src.length)
          throw new JslSyntaxError('Unclosed string', pos);
        if (src[i] === '"')
          break;
        if (src.startsWith('\\!', i) && i + 2 < src.length) {
          text += ESCAPES[src[i + 2]] ?? src[i + 2];
          i += 3;
          continue;
        }
        text += src[i++];
      }
      i++;
      if (src[i] === 'n' && !/[\w]/.test(src[i + 1] ?? '')) {
        tokens.push({kind: 'nameLiteral', text, pos});
        i++;
      } else
        tokens.push({kind: 'string', text, pos});
      continue;
    }
    // a lone `.` is JSL's missing value
    if (c === '.' && !/[0-9]/.test(src[i + 1] ?? '')) {
      tokens.push({kind: 'number', text: '.', pos});
      i++;
      continue;
    }
    if (/[0-9]/.test(c) || (c === '.' && /[0-9]/.test(src[i + 1] ?? ''))) {
      const m = /^(\d*\.?\d+(?:[eE][+-]?\d+)?|\d+\.)/.exec(src.slice(i))!;
      // dates such as 01Jan2020 continue as a name
      if (/[A-Za-z_]/.test(src[i + m[0].length] ?? '')) {
        const n = /^[\w.]+/.exec(src.slice(i))!;
        tokens.push({kind: 'name', text: n[0], pos});
        i += n[0].length;
      } else {
        tokens.push({kind: 'number', text: m[0], pos});
        i += m[0].length;
      }
      continue;
    }
    if (/[A-Za-z_]/.test(c)) {
      // a name may contain spaces between words: `Fit Y by X`
      const m = /^[A-Za-z_][\w]*(?:[ \t]+[A-Za-z_\d][\w]*)*/.exec(src.slice(i))!;
      tokens.push({kind: 'name', text: m[0], pos});
      i += m[0].length;
      continue;
    }
    const p = PUNCT.find((s) => src.startsWith(s, i));
    if (!p)
      throw new JslSyntaxError(`Unexpected character "${c}"`, pos);
    tokens.push({kind: 'punct', text: p, pos});
    i += p.length;
  }
  tokens.push({kind: 'end', text: '', pos: src.length});
  return tokens;
}

// Binary operators by precedence, loosest first. `;` glues statements; `<<` sends a message.
const BINARY_LEVELS: string[][] = [
  [';'],
  ['=', '+=', '-=', '*=', '/='],
  ['||', '|'],
  ['&&', '&'],
  ['==', '!=', '<', '>', '<=', '>='],
  ['<<'],
  ['+', '-'],
  ['*', '/'],
  ['^'],
];
const MAX_DEPTH = 200;

/** Parses JSL source into a tree. Throws a {@link JslSyntaxError} for text it can't read. */
export function parseJsl(src: string): JslNode {
  const tokens = tokenize(src);
  let k = 0;
  let depth = 0;

  const peek = (): Token => tokens[k];
  const isPunct = (text: string): boolean => tokens[k].kind === 'punct' && tokens[k].text === text;
  const expectPunct = (text: string): void => {
    if (!isPunct(text))
      throw new JslSyntaxError(`Expected "${text}"`, tokens[k].pos);
    k++;
  };

  function args(close: string): JslNode[] {
    if (++depth > MAX_DEPTH)
      throw new JslSyntaxError('Script is nested too deeply', peek().pos);
    const out: JslNode[] = [];
    while (!isPunct(close)) {
      if (peek().kind === 'end')
        throw new JslSyntaxError(`Expected "${close}"`, peek().pos);
      // empty arguments and trailing commas are allowed
      if (isPunct(',')) {
        k++;
        continue;
      }
      out.push(binary(1));
      if (!isPunct(close))
        expectPunct(',');
    }
    k++;
    depth--;
    return out;
  }

  function binary(level: number): JslNode {
    if (level >= BINARY_LEVELS.length)
      return unary();
    let left = binary(level + 1);
    while (peek().kind === 'punct' && BINARY_LEVELS[level].includes(peek().text)) {
      const op = tokens[k++].text;
      // a trailing `;` ends the script
      if (op === ';' && (peek().kind === 'end' || isPunct(')') || isPunct(',')))
        break;
      left = {type: 'binary', op, left, right: binary(level + 1)};
    }
    return left;
  }

  function unary(): JslNode {
    // a leading `<<` sends to the object before it, as in `Fit Group(..., <<{Profiler(1)})`
    if (isPunct('-') || isPunct('!') || isPunct('+') || isPunct('<<')) {
      if (++depth > MAX_DEPTH)
        throw new JslSyntaxError('Script is nested too deeply', peek().pos);
      const op = tokens[k++].text;
      const operand = unary();
      depth--;
      return {type: 'unary', op, operand};
    }
    return postfix(primary());
  }

  function postfix(node: JslNode): JslNode {
    for (;;) {
      if (isPunct('[')) {
        k++;
        node = {type: 'index', target: node, args: args(']')};
      } else if (isPunct('(') && isCall(node, 'Name') && node.args.length === 1 && node.args[0].type === 'string') {
        // `Name("1a fit")(...)` calls the name `1a fit`, as JMP writes a script named so
        k++;
        node = {type: 'call', name: node.args[0].value, args: args(')')};
      } else
        return node;
    }
  }

  function primary(): JslNode {
    const t = tokens[k];
    switch (t.kind) {
    case 'number':
      k++;
      return {type: 'number', value: t.text === '.' ? NaN : Number(t.text)};
    case 'string':
      k++;
      // `"Value Ordering"(...)`: JMP writes a property it has no name for as a string called
      if (isPunct('(')) {
        k++;
        return {type: 'call', name: t.text, args: args(')')};
      }
      return {type: 'string', value: t.text};
    case 'nameLiteral':
    case 'name':
      k++;
      if (isPunct('(')) {
        k++;
        return {type: 'call', name: t.text, args: args(')')};
      }
      return {type: 'name', name: t.text};
    case 'punct':
      if (t.text === '(') {
        k++;
        if (++depth > MAX_DEPTH)
          throw new JslSyntaxError('Script is nested too deeply', t.pos);
        const inner = isPunct(')') ? {type: 'list', items: []} as JslNode : binary(0);
        expectPunct(')');
        depth--;
        return inner;
      }
      if (t.text === '{') {
        k++;
        return {type: 'list', items: args('}')};
      }
      if (t.text === '[') {
        k++;
        return {type: 'list', items: args(']')};
      }
      if (t.text === ':' || t.text === '::') {
        k++;
        return column();
      }
      throw new JslSyntaxError(`Unexpected "${t.text}"`, t.pos);
    default:
      throw new JslSyntaxError('Unexpected end of script', t.pos);
    }
  }

  /** After `:`: a column reference such as `:Age`, `:"Age at start"n` or `:Name("Age at start")`. */
  function column(): JslNode {
    const t = tokens[k];
    if (t.kind === 'nameLiteral') {
      k++;
      return {type: 'column', name: t.text};
    }
    if (t.kind === 'name') {
      k++;
      if (jslKey(t.text) === 'name' && isPunct('(')) {
        k++;
        const inner = args(')');
        if (inner.length === 1 && inner[0].type === 'string')
          return {type: 'column', name: inner[0].value};
        throw new JslSyntaxError('Expected a column name in Name()', t.pos);
      }
      return {type: 'column', name: t.text};
    }
    throw new JslSyntaxError('Expected a column name after ":"', t.pos);
  }

  const root = peek().kind === 'end' ? {type: 'list', items: []} as JslNode : binary(0);
  if (peek().kind !== 'end')
    throw new JslSyntaxError(`Unexpected "${peek().text}"`, peek().pos);
  return root;
}

/** Statements of a script, flattening `a; b; c`. */
export function statements(node: JslNode): JslNode[] {
  return node.type === 'binary' && node.op === ';' ? [...statements(node.left), ...statements(node.right)] : [node];
}

/** True when `node` is a call to `name` (compared the way JSL compares names). */
export function isCall(node: JslNode | undefined, name: string): node is Extract<JslNode, {type: 'call'}> {
  return node?.type === 'call' && jslKey(node.name) === jslKey(name);
}

/** The first argument of `node` that is a call to `name`. */
export function findArg(node: Extract<JslNode, {type: 'call'}>, name: string):
  Extract<JslNode, {type: 'call'}> | undefined {
  return node.args.find((a) => isCall(a, name)) as Extract<JslNode, {type: 'call'}> | undefined;
}

/** Every call to `name` anywhere under `node`, outermost first. */
export function findAll(node: JslNode, name: string): Extract<JslNode, {type: 'call'}>[] {
  const out: Extract<JslNode, {type: 'call'}>[] = [];
  const visit = (n: JslNode): void => {
    if (isCall(n, name))
      out.push(n);
    for (const child of children(n))
      visit(child);
  };
  visit(node);
  return out;
}

function children(n: JslNode): JslNode[] {
  switch (n.type) {
  case 'call': return n.args;
  case 'list': return n.items;
  case 'index': return [n.target, ...n.args];
  case 'unary': return [n.operand];
  case 'binary': return [n.left, n.right];
  default: return [];
  }
}

/** Column names given to a role, as in `Y(:Weight, :Height)` or `Column("Age")`. */
export function columnNames(node: Extract<JslNode, {type: 'call'}> | undefined): string[] {
  if (!node)
    return [];
  const out: string[] = [];
  for (const a of node.args) {
    if (a.type === 'column')
      out.push(a.name);
    else if (a.type === 'string')
      out.push(a.value);
    else if (a.type === 'list')
      out.push(...columnNames({type: 'call', name: '', args: a.items}));
    else if (isCall(a, 'Column') && a.args.length === 1 && a.args[0].type === 'string')
      out.push(a.args[0].value);
  }
  return out;
}
