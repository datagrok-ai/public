import {FilterError} from './model.js';
import type {DomainCondition, DomainConditionNode, DomainConditionTree, FilterProblem} from './model.js';
import {property as propertyOf} from './schema.js';
import type {FilterProperty, FilterSchema} from './schema.js';
import {operators} from './operators.js';
import type {FilterOperator} from './operators.js';
import {markSpan, resolveSpan, spanOf} from '../span.js';

export interface FilterToken {
  type: 'name' | 'tag' | 'string' | 'number' | 'symbol' | 'punct';
  text: string;
  start: number;
  end: number;
}

/** What the grammar expects at the caret — the query input's completion driver. */
export interface FilterCompletion {
  expect: 'property' | 'operator' | 'value' | 'connector';
  /** The text of the token under the caret up to the caret (a string's without its quote). */
  prefix: string;
  replace: {start: number, end: number};
  property?: FilterProperty;
  operator?: FilterOperator;
}

const WS = /[\t\n\v\f\r \u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000\ufeff]/;
const LETTER = /[A-Za-z]/;
const WORD = /[A-Za-z0-9_]/;
const DIGIT = /[0-9]/;
const HEX4 = /^[0-9A-Fa-f]{4}/;
const NUMBER = /^-?(?:0|\d+)(?:\.\d+)?(?:[eE][-+]?\d+)?/;
const INTEGER = /^-?\d+$/;
const ESCAPES: Record<string, string> = {
  '\\': '\\', '/': '/', '"': '"', '\'': '\'', 'b': '\b', 'f': '\f', 'n': '\n', 'r': '\r', 't': '\t',
};
const OPERATORS = ['!matches', '!like', 'matches', 'starts', 'like', 'ends', '!~', '~=', '!=', '>=', '<=', '=', '~',
  '>', '<'];
const SYMBOLS = ['!matches', '!like', '!in', '!~', '~=', '!=', '>=', '<=', '&&', '||', '=', '~', '>', '<', '&', '|'];
const NUMBER_SPANS: [string, number][] = [['kb', 1024], ['KB', 1024], ['MB', 1024 ** 2], ['GB', 1024 ** 3],
  ['TB', 1024 ** 4], ['k', 1000], ['K', 1000], ['M', 1e6]];
const FLIPPED: Record<string, string> = {'<': '>', '>': '<', '<=': '>=', '>=': '<='};
const NEGATED: Record<string, string> = {
  '=': '!=', '!=': '=', '>': '<=', '<=': '>', '<': '>=', '>=': '<', 'like': 'not like', 'not like': 'like',
  '~*': '!~*', '!~*': '~*', 'is': 'is not', 'is not': 'is',
};
const LITERALS = new Set(['true', 'false', 'null', 'now', '@current']);
const WORD_OPERATORS = new Set(['like', 'starts', 'ends', 'matches', 'fuzzy', 'in', 'between']);
const SYMBOL_OPERATORS = new Set(['=', '!=', '>', '>=', '<', '<=', '~', '~=', '!~', '!like', '!matches', '!in']);
const SPELLINGS: Record<string, string> = {'~': 'matches', '!~': '!matches', '~=': 'fuzzy', '!in': 'not in'};
const CONNECTORS = new Set(['and', 'or', '&&', '&', '||', '|']);

const node = (property: string, operator: string, value: unknown): DomainCondition => ({property, operator, value});
export const isCondition = (n: DomainConditionNode): n is DomainCondition =>
  typeof n === 'object' && !Array.isArray(n);
const dartDouble = (n: number): string => Number.isInteger(n) && Math.abs(n) < 1e21 ? `${n}.0` : String(n);

/** `[{fuzzy}, 'or', {like}]` on one property, or null. */
export function fuzzyPair(n: DomainConditionNode): {fuzzy: DomainCondition, like: DomainCondition} | null {
  if (!Array.isArray(n) || n.length !== 3 || n[1] !== 'or' || !isCondition(n[0]) || !isCondition(n[2]))
    return null;
  const [fuzzy, , like] = n as [DomainCondition, 'or', DomainCondition];
  return fuzzy.operator === 'fuzzy' && like.operator === 'like' && fuzzy.property === like.property ?
    {fuzzy, like} : null;
}

/** De Morgan pushdown: connectors swap, operators flip, lists keep their list, a
 * fuzzy pair becomes the single `not like`. */
export function negate(tree: DomainConditionTree): DomainConditionTree {
  return tree.map((n): DomainConditionNode => {
    if (n === 'and')
      return 'or';
    if (n === 'or')
      return 'and';
    if (Array.isArray(n)) {
      const pair = fuzzyPair(n);
      return pair ? node(pair.like.property, 'not like', pair.like.value) : negate(n);
    }
    const operator = NEGATED[n.operator];
    if (!operator)
      throw new FilterError(`Cannot negate "${n.operator}"`);
    return {...n, operator};
  });
}

/** Ordered-choice recursive descent over the raw text, one method per Dart production; every
 * token skips whitespace on both sides (`token()`), failures rewind and record the farthest
 * expectation for the syntax problem. */
class Reader {
  pos = 0;
  farthest = 0;
  expected: string[] = [];
  /** Dart's spelling of the number `value()` last read: `5.0` stays a double, `5k` is `5000.0`. */
  numberText = '';

  constructor(readonly text: string, readonly now: Date) {}

  /** Dart's `'$value'` — what `like` and `fuzzy` interpolate into their patterns. */
  interpolate(value: unknown): string {
    return value instanceof Date ? value.toISOString().replace('T', ' ') :
      typeof value === 'number' ? this.numberText : String(value);
  }

  fail(what: string, at: number = this.pos): undefined {
    if (at > this.farthest) {
      this.farthest = at;
      this.expected = [];
    }
    if (at === this.farthest && !this.expected.includes(what))
      this.expected.push(what);
    return undefined;
  }

  ws(): void {
    while (this.pos < this.text.length && WS.test(this.text[this.pos]))
      this.pos++;
  }

  ch(offset: number = 0): string {
    return this.text[this.pos + offset] ?? '';
  }

  literal(s: string): boolean {
    if (!this.text.startsWith(s, this.pos))
      return false;
    this.pos += s.length;
    return true;
  }

  /** `token(string(s))`: whitespace-trimmed on both sides. */
  token(s: string): boolean {
    const save = this.pos;
    this.ws();
    if (this.literal(s)) {
      this.ws();
      return true;
    }
    this.pos = save;
    return false;
  }

  /** `token(string(s) & word().not())`. */
  keyword(s: string): boolean {
    const save = this.pos;
    this.ws();
    if (this.literal(s) && !WORD.test(this.ch())) {
      this.ws();
      return true;
    }
    this.pos = save;
    return false;
  }

  start(): DomainConditionTree | undefined {
    this.ws();
    if (this.pos === this.text.length)
      return [];
    const tree = this.conditions();
    if (tree === undefined)
      return undefined;
    this.ws();
    return this.pos === this.text.length ? tree : this.fail('end of input');
  }

  /** `term ((connector | ε) term)*`; a list term is unwrapped when it holds one element — the
   * Dart `conditions()` map — except a parenthesized one, which always nests. */
  conditions(): DomainConditionTree | undefined {
    const first = this.term();
    if (first === undefined)
      return undefined;
    const list: DomainConditionTree = [first];
    for (;;) {
      const save = this.pos;
      const connector = this.connector();
      const next = this.term();
      if (next === undefined) {
        this.pos = save;
        return list;
      }
      list.push(connector, next);
    }
  }

  connector(): 'and' | 'or' {
    if (this.keyword('or') || this.token('||') || this.token('|'))
      return 'or';
    this.keyword('and') || this.token('&&') || this.token('&');
    return 'and';
  }

  term(): DomainConditionNode | undefined {
    this.ws();
    const start = this.pos;
    const condition = this.condition();
    if (condition !== undefined)
      return condition;
    this.pos = start;
    if (this.keyword('not') && this.token('(')) {
      const inner = this.conditions();
      if (inner !== undefined && this.token(')')) {
        const negated = negate(inner);
        return negated.length === 1 ? negated[0] : negated;
      }
      if (inner !== undefined)
        this.fail('")"');
    }
    this.pos = start;
    if (this.token('(')) {
      const inner = this.conditions();
      if (inner !== undefined && this.token(')'))
        return inner;
      if (inner !== undefined)
        this.fail('")"');
    }
    this.pos = start;
    return this.fail('a condition');
  }

  condition(): DomainConditionNode | undefined {
    const start = this.pos;
    const property = this.property();
    if (property !== undefined) {
      const after = this.pos;
      for (const alternative of [this.fuzzyCondition, this.operatorCondition, this.inCondition,
        this.betweenCondition]) {
        const result = alternative.call(this, property);
        if (result !== undefined)
          return result;
        this.pos = after;
      }
    }
    this.pos = start;
    if (this.keyword('not')) {
      const negated = this.property();
      if (negated !== undefined)
        return node(negated, '=', false);
    }
    this.pos = start;
    this.ws();
    if (this.literal('#')) {
      const from = this.pos;
      while (WORD.test(this.ch()))
        this.pos++;
      const tag = this.text.slice(from, this.pos);
      this.ws();
      return node('entityTags.tag', '=', tag);
    }
    this.pos = start;
    return undefined;
  }

  /** `(name | [text]) ('.' …)*`, flattened — brackets and their escapes travel verbatim. */
  property(): string | undefined {
    this.ws();
    const start = this.pos;
    if (!this.nameSegment()) {
      this.pos = start;
      return undefined;
    }
    while (this.ch() === '.') {
      const save = this.pos;
      this.pos++;
      if (!this.nameSegment()) {
        this.pos = save;
        break;
      }
    }
    const name = this.text.slice(start, this.pos);
    this.ws();
    return name;
  }

  nameSegment(): boolean {
    if (LETTER.test(this.ch())) {
      this.pos++;
      while (WORD.test(this.ch()))
        this.pos++;
      return true;
    }
    if (this.ch() !== '[')
      return false;
    this.pos++;
    while (this.pos < this.text.length) {
      const c = this.ch();
      if (c === ']') {
        this.pos++;
        return true;
      }
      if (c !== '\\')
        this.pos++;
      else if (this.ch(1) === ']' || this.ch(1) in ESCAPES)
        this.pos += 2;
      else if (this.ch(1) === 'u' && HEX4.test(this.text.slice(this.pos + 2, this.pos + 6)))
        this.pos += 6;
      else
        return false;
    }
    return false;
  }

  operator(): string | undefined {
    this.ws();
    for (const op of OPERATORS) {
      if (this.literal(op)) {
        this.ws();
        return op;
      }
    }
    return undefined;
  }

  fuzzyCondition(property: string): DomainConditionNode | undefined {
    if (!this.keyword('fuzzy'))
      return undefined;
    let threshold: number | null = null;
    if (this.ch() === '(') {
      const save = this.pos;
      this.pos++;
      this.ws();
      const n = this.number();
      this.ws();
      if (n !== undefined && this.literal(')'))
        threshold = n;
      else
        this.pos = save;
    }
    const value = this.value();
    if (value === undefined)
      return undefined;
    return [{property, operator: 'fuzzy', threshold, value}, 'or',
      node(property, 'like', `%${this.interpolate(value)}%`)];
  }

  operatorCondition(property: string): DomainConditionNode | undefined {
    const operator = this.operator();
    if (operator === undefined)
      return this.fail('an operator');
    const value = this.value();
    return value === undefined ? undefined : this.normalize(property, operator, value);
  }

  /** The Dart `condition()` map: span flip, `like`/`starts`/`ends`/`!like` shapes, regex and
   * fuzzy spellings. */
  normalize(property: string, operator: string, value: unknown): DomainConditionNode {
    const span = value instanceof Date ? spanOf(value) : undefined;
    // Dart flips a span only when its date lands after now + 1 s, so `0h` and `-0h` stay put
    if (span !== undefined && (value as Date).getTime() > this.now.getTime() + 1000) {
      value = markSpan(resolveSpan(`-${span}`, this.now), `-${span}`);
      operator = FLIPPED[operator] ?? operator;
    }
    const text = this.interpolate(value);
    switch (operator) {
      case 'like': return node(property, 'like', `%${text}%`);
      case 'starts': return node(property, 'like', `${text}%`);
      case 'ends': return node(property, 'like', `%${text}`);
      case '!like': return node(property, 'not like', `%${text}%`);
      case 'matches': case '~': return node(property, '~*', value);
      case '!matches': case '!~': return node(property, '!~*', value);
      case '~=':
        return [{property, operator: 'fuzzy', threshold: null, value}, 'or', node(property, 'like', `%${text}%`)];
      default: return node(property, operator, value);
    }
  }

  inCondition(property: string): DomainConditionNode | undefined {
    const operator = this.keyword('in') ? '=' : this.keyword('not in') || this.keyword('!in') ? '!=' : undefined;
    if (operator === undefined)
      return undefined;
    if (!this.token('('))
      return this.fail('"("');
    const values: unknown[] = [];
    for (;;) {
      const value = this.value();
      if (value === undefined)
        return undefined;
      values.push(value);
      if (this.token(')'))
        return node(property, operator, values);
      if (!this.token(','))
        return this.fail('"," or ")"');
    }
  }

  betweenCondition(property: string): DomainConditionNode | undefined {
    if (!this.keyword('between'))
      return undefined;
    const low = this.value();
    if (low === undefined)
      return undefined;
    if (!this.keyword('and'))
      return this.fail('"and"');
    const high = this.value();
    if (high === undefined)
      return undefined;
    return [node(property, '>=', low), 'and', node(property, '<=', high)];
  }

  /** `true | false | null | @current | timeSpan | numberSpan | number | string`, in that order. */
  value(): unknown {
    this.ws();
    const start = this.pos;
    const done = (v: unknown): unknown => {
      this.ws();
      return v;
    };
    if (this.literal('true'))
      return done(true);
    if (this.literal('false'))
      return done(false);
    if (this.literal('null'))
      return done(null);
    if (this.literal('@current'))
      return done('@current');
    const m = NUMBER.exec(this.text.slice(start));
    if (m) {
      const unit = this.text[start + m[0].length] ?? '';
      if ('hwdmy'.includes(unit) && unit !== '') {
        if (!INTEGER.test(m[0]))
          return this.fail('a whole number of h, d, w, m or y', start);
        this.pos = start + m[0].length + 1;
        const span = `${m[0]}${unit}`;
        return done(markSpan(resolveSpan(span, this.now), span));
      }
    }
    if (this.literal('now'))
      return done(markSpan(new Date(this.now.getTime()), 'now'));
    if (m) {
      this.pos = start + m[0].length;
      for (const [suffix, factor] of NUMBER_SPANS) {
        if (this.literal(suffix)) {
          const n = Number(m[0]) * factor;
          this.numberText = dartDouble(n);
          return done(n);
        }
      }
      const n = Number(m[0]);
      this.numberText = m[0].includes('.') ? dartDouble(n) : String(n);
      // `-0` → 0: Dart reads `-0` as the integer 0, and deep equality tells the two apart
      return done(n === 0 ? 0 : n);
    }
    const s = this.string();
    return s === undefined ? this.fail('a value', start) : done(s);
  }

  number(): number | undefined {
    const m = NUMBER.exec(this.text.slice(this.pos));
    if (!m)
      return undefined;
    this.pos += m[0].length;
    return Number(m[0]);
  }

  /** The quote that opens is the one that closes; `\` escapes and `\uXXXX`. */
  string(): string | undefined {
    const quote = this.ch();
    if (quote !== '"' && quote !== '\'')
      return undefined;
    this.pos++;
    let out = '';
    for (;;) {
      if (this.pos >= this.text.length)
        return this.fail('a closing quote');
      const c = this.ch();
      if (c === quote) {
        this.pos++;
        return out;
      }
      if (c !== '\\') {
        out += c;
        this.pos++;
      } else if (this.ch(1) in ESCAPES) {
        out += ESCAPES[this.ch(1)];
        this.pos += 2;
      } else if (this.ch(1) === 'u' && HEX4.test(this.text.slice(this.pos + 2, this.pos + 6))) {
        out += String.fromCharCode(parseInt(this.text.slice(this.pos + 2, this.pos + 6), 16));
        this.pos += 6;
      } else
        return this.fail('an escape');
    }
  }
}

/** The Dart tree for `text` (`errors` empty), or the syntax problems with positions. Blank
 * text is the empty tree (the server's reading). */
export function parseTree(text: string, options?: {now?: Date}):
  {tree: DomainConditionTree, errors: FilterProblem[]} {
  const reader = new Reader(text, options?.now ?? new Date());
  const tree = reader.start();
  if (tree !== undefined)
    return {tree, errors: []};
  const start = reader.farthest;
  const token = tokenize(text).find((t) => t.start <= start && start < t.end);
  const end = token ? token.end : Math.min(start + 1, text.length);
  return {tree: [], errors: [{
    nodeId: null, code: 'syntax', message: `Expected ${reader.expected.join(' or ')}`, position: {start, end},
  }]};
}

function tokenize(text: string): FilterToken[] {
  const tokens: FilterToken[] = [];
  let pos = 0;
  const push = (type: FilterToken['type'], start: number) =>
    tokens.push({type, text: text.slice(start, pos), start, end: pos});
  const skipEscaped = () => pos += text[pos] === '\\' ? 2 : 1;
  while (pos < text.length) {
    const c = text[pos];
    const start = pos;
    if (WS.test(c))
      pos++;
    else if (c === '"' || c === '\'') {
      pos++;
      while (pos < text.length && text[pos] !== c)
        skipEscaped();
      pos = Math.min(pos + 1, text.length);
      push('string', start);
    } else if (c === '#' || c === '@') {
      pos++;
      while (WORD.test(text[pos] ?? ''))
        pos++;
      push(c === '#' ? 'tag' : 'name', start);
    } else if (LETTER.test(c) || c === '_' || c === '[') {
      for (;;) {
        if (text[pos] === '[') {
          pos++;
          while (pos < text.length && text[pos] !== ']')
            skipEscaped();
          pos = Math.min(pos + 1, text.length);
        } else {
          while (WORD.test(text[pos] ?? ''))
            pos++;
        }
        if (text[pos] === '.' && (WORD.test(text[pos + 1] ?? '') || text[pos + 1] === '['))
          pos++;
        else
          break;
      }
      push('name', start);
    } else if (DIGIT.test(c) || (c === '-' && DIGIT.test(text[pos + 1] ?? ''))) {
      pos += NUMBER.exec(text.slice(pos))![0].length;
      while (LETTER.test(text[pos] ?? ''))
        pos++;
      push('number', start);
    } else if (c === '(' || c === ')' || c === ',') {
      pos++;
      push('punct', start);
    } else {
      pos += SYMBOLS.find((s) => text.startsWith(s, pos))?.length ?? 1;
      push('symbol', start);
    }
  }
  return tokens;
}

export function completionContext(text: string, caret: number, schema: FilterSchema): FilterCompletion {
  const tokens = tokenize(text);
  const closed = (t: FilterToken) => t.type === 'string' && t.end - t.start > 1 && t.text.endsWith(t.text[0]);
  const current = tokens.find((t) => t.type !== 'punct' && t.start < caret && caret <= t.end &&
    !(caret === t.end && closed(t)));
  let expect: FilterCompletion['expect'] = 'property';
  let property: FilterProperty | undefined;
  let operatorId: string | undefined;
  let list = false;
  let between = 0;
  let threshold = false;
  let pendingNot = false;
  const isValue = (t: FilterToken) => t.type === 'string' || t.type === 'number' ||
    (t.type === 'name' && LITERALS.has(t.text));
  const startCondition = (t: FilterToken) => {
    property = propertyOf(schema, t.text) ?? undefined;
    operatorId = undefined;
  };
  for (const t of tokens) {
    if (t === current || t.end > caret)
      break;
    const word = t.type === 'name' ? t.text : '';
    if (expect === 'property') {
      if (t.type === 'tag')
        expect = 'connector';
      else if (word === 'not')
        pendingNot = true;
      else if (t.type === 'name' && !LITERALS.has(word)) {
        startCondition(t);
        expect = pendingNot ? 'connector' : 'operator';
        pendingNot = false;
      } else
        pendingNot = false;
    } else if (expect === 'operator') {
      if (word === 'not') {
        pendingNot = true;
        continue;
      }
      const spelling = pendingNot && word === 'in' ? 'not in' : t.text;
      pendingNot = false;
      if ((t.type === 'symbol' && SYMBOL_OPERATORS.has(t.text)) || (t.type === 'name' && WORD_OPERATORS.has(word))) {
        operatorId = SPELLINGS[spelling] ?? spelling;
        list = operatorId === 'in' || operatorId === 'not in';
        between = operatorId === 'between' ? 1 : 0;
        expect = 'value';
      }
    } else if (expect === 'value') {
      if (t.type === 'punct' && t.text === '(') {
        if (operatorId === 'fuzzy')
          threshold = true;
      } else if (t.type === 'punct' && t.text === ')') {
        if (threshold)
          threshold = false;
        else if (list) {
          list = false;
          expect = 'connector';
        }
      } else if (between === 1 && word === 'and')
        between = 2;
      else if (isValue(t) && !list && !threshold && between !== 1) {
        between = 0;
        expect = 'connector';
      }
    } else if (((t.type === 'symbol' || t.type === 'name') && CONNECTORS.has(t.text)) || word === 'not') {
      expect = 'property';
      pendingNot = word === 'not';
      property = undefined;
      operatorId = undefined;
    } else if (t.type === 'name' && !LITERALS.has(word)) {
      startCondition(t);
      expect = 'operator';
    } else if (t.type === 'tag')
      expect = 'connector';
  }
  const prefix = !current ? '' : text.slice(current.start + (current.type === 'string' ? 1 : 0), caret);
  const replace = current ? {start: current.start, end: current.end} : {start: caret, end: caret};
  const result: FilterCompletion = {expect, prefix, replace};
  if (property)
    result.property = property;
  const operator = operatorId === undefined ? undefined :
    (property && operators.get(operatorId, property)) ?? operators.get(operatorId);
  if (operator)
    result.operator = operator;
  return result;
}
