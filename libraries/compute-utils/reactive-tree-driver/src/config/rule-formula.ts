import {Grammars, IToken, Parser} from 'ebnf';
import {RestrictionType} from '../data/common-types';
import {driverOpNames, isLiteral} from '../runtime/rule-expressions';
import {
  PipelineCheckConfiguration, PipelineRuleConfiguration, RuleEffect, RuleExpr, RuleLogic, RuleSource, SourceKind,
} from './PipelineConfiguration';

/* eslint-disable max-len */
const formulaGrammar = `
Expr     ::= WS* Value WS* {fragment=true}
Value    ::= Call | Array | String | Number | Ref {fragment=true}
Call     ::= Name WS* '(' WS* (Arg (WS* ',' WS* Arg)*)? WS* ')'
Arg      ::= NamedArg | Expr {fragment=true}
NamedArg ::= Name WS* ':' WS* Expr
Array    ::= '[' WS* (Expr (WS* ',' WS* Expr)*)? WS* ']'
Ref      ::= RefRoot ('.' Segment)*
RefRoot  ::= '$'? Name
Segment  ::= Name | Digits
String   ::= '"' Char* '"' | #x60 Raw? #x60
Char     ::= [^"#x5C] | #x5C [#x0-#xFFFF] {fragment=true}
Raw      ::= [^#x60]+ {fragment=true}
Number   ::= '-'? (Digits ('.' Digits?)? | '.' Digits) (('e' | 'E') ('+' | '-')? Digits)?
Digits   ::= [0-9]+
Name     ::= [_a-zA-Z][a-zA-Z_0-9]*
WS       ::= [#x20#x09#x0A#x0D]
`;
/* eslint-enable max-len */

const formulaParser = new Parser(Grammars.Custom.getRules(formulaGrammar));

// names come from user text, so lookups must not see Object.prototype (toString, constructor, ...)
const has = (table: object, key: string) => Object.hasOwn(table, key);

// JSON Logic ops whose names are not identifiers
const renamed: Record<string, string> = {
  not: '!', bool: '!!', eq: '==', ne: '!=', same: '===', notSame: '!==', gt: '>', gte: '>=', lt: '<', lte: '<=',
  add: '+', sub: '-', mul: '*', div: '/', mod: '%',
};

/** Ops a formula calls by their own name: JSON Logic's and the driver's. */
export const formulaOps = new Set([
  'if', 'and', 'or', 'in', 'cat', 'substr', 'log', 'min', 'max', 'merge',
  'map', 'filter', 'reduce', 'all', 'some', 'none', 'var', 'missing', 'missing_some',
  ...driverOpNames,
]);

// the element of these ops' second argument is `$it`; `reduce` reads `current` and `accumulator` instead
const elementOps = new Set(['map', 'filter', 'all', 'some', 'none']);
const ELEMENT = '$it';

type EffectName = RuleEffect['effect'];

// the effect object behind a name; found by membership because hide and show share one shape
type EffectOf<E extends EffectName> = RuleEffect extends infer R ? R extends {effect: infer N} ?
  E extends N ? R : never : never : never;
// what a call has to cover: every field except the name and the condition every call takes
type EffectFields<E extends EffectName> = Exclude<keyof EffectOf<E>, 'effect' | 'when'>;
type EffectCall<E extends EffectName> = {positional: readonly EffectFields<E>[], options: readonly EffectFields<E>[]};

// positional parameters and named options of each effect call;
// `meta` takes any other named argument as a meta key, which is what its `meta` option stands for
const effectCalls = {
  hide: {positional: ['targets'], options: []},
  show: {positional: ['targets'], options: []},
  items: {positional: ['targets', 'items'], options: []},
  meta: {positional: ['targets'], options: ['meta']},
  error: {positional: ['targets', 'message'], options: []},
  warning: {positional: ['targets', 'message'], options: []},
  notification: {positional: ['targets', 'message'], options: []},
  verdicts: {positional: ['targets', 'source'], options: []},
  set: {positional: ['targets', 'value'], options: ['restriction']},
  clear: {positional: ['targets'], options: ['restriction']},
  assign: {positional: ['values'], options: ['targets', 'restriction', 'ignoreCase']},
} as const satisfies {[E in EffectName]: EffectCall<E>};

// a field of an effect type that its call does not list fails the build, naming it as "effect.field"
type ListedFields<E extends EffectName> =
  (typeof effectCalls)[E]['positional'][number] | (typeof effectCalls)[E]['options'][number];
type UnlistedFields = {[E in EffectName]: Exclude<EffectFields<E>, ListedFields<E>> extends infer M ?
  [M] extends [never] ? never : `${E}.${M & string}` : never}[EffectName];
type NoneUnlisted<T extends never> = T;
type _EveryEffectFieldHasACall = NoneUnlisted<UnlistedFields>;

// `js` holds a function and a dataframe `table` a live object, so both stay objects
const sourceCalls: Record<Exclude<SourceKind, 'js'>, {positional: string[]}> = {
  file: {positional: ['path']},
  table: {positional: ['csv']},
  func: {positional: ['name']},
  query: {positional: ['connection', 'sql']},
  validators: {positional: ['input']},
  choices: {positional: ['input']},
};

const restrictions: Record<RestrictionType, true> = {disabled: true, restricted: true, info: true, none: true};

class FormulaError extends Error {
  constructor(message: string, node: {start?: number} | undefined, text: string) {
    super(`${message} at column ${(node?.start ?? 0) + 1}: ${text}`);
  }
}

function parse(text: string): IToken {
  const ast = formulaParser.getAST(text, 'Expr');
  if (!ast || ast.rest || ast.end !== text.length)
    throw new FormulaError('Syntax error', {start: ast?.end ?? 0}, text);
  return ast.type === 'Expr' ? ast.children[0] : ast;
}

const escapes: Record<string, string> = {'n': '\n', 't': '\t', '"': '"', '\\': '\\'};

function unquote(node: IToken, text: string): string {
  const raw = node.text;
  if (raw[0] === '`')
    return raw.slice(1, -1);
  return raw.slice(1, -1).replace(/\\(.)/gs, (escape, char, offset) => {
    if (!has(escapes, char)) {
      throw new FormulaError(`Unknown escape ${escape}, write patterns and other text with backslashes in backticks`,
        {start: node.start + 1 + offset}, text);
    }
    return escapes[char];
  });
}

function callParts(node: IToken, text: string) {
  const [name, ...args] = node.children;
  const positional: IToken[] = [];
  const named: Record<string, IToken> = {};
  for (const arg of args) {
    if (arg.type === 'NamedArg') {
      const key = arg.children[0].text;
      if (has(named, key))
        throw new FormulaError(`Duplicate argument ${key}`, arg, text);
      named[key] = arg.children[1];
    } else {
      if (Object.keys(named).length)
        throw new FormulaError('Positional argument after a named one', arg, text);
      positional.push(arg);
    }
  }
  return {name: name.text, positional, named};
}

function isConstant(value: any): boolean {
  if (Array.isArray(value))
    return value.every(isConstant);
  return value == null || typeof value !== 'object' || isLiteral(value);
}

function unliteral(value: any): any {
  if (Array.isArray(value))
    return value.map(unliteral);
  return isLiteral(value) ? value.literal : value;
}

// missing, missing_some and var take key names, so a bare alias there is its path, not a lookup
function pathArg(node: IToken | undefined, text: string, op: string): string {
  if (node?.type === 'Ref' && (node.text === ELEMENT || node.text.startsWith(`${ELEMENT}.`)))
    throw new FormulaError(`${op} takes aliases, read the element as ${ELEMENT} and its fields by name`, node, text);
  if (node?.type === 'Ref')
    return node.text;
  if (node?.type === 'String')
    return unquote(node, text);
  throw new FormulaError(`${op} takes aliases`, node, text);
}

function expression(node: IToken, text: string, inElement = false): any {
  switch (node.type) {
  case 'Number': return Number(node.text);
  case 'String': return unquote(node, text);
  case 'Array': return node.children.map((child) => expression(child, text, inElement));
  case 'Ref': {
    const path = node.text;
    if (path === 'true' || path === 'false')
      return path === 'true';
    if (path === 'null')
      return null;
    if (path.startsWith(`${ELEMENT}.`))
      throw new FormulaError(`Read the element's fields by name, ${path.slice(ELEMENT.length + 1)}`, node, text);
    if (path === ELEMENT) {
      if (!inElement)
        throw new FormulaError(`${ELEMENT} is the element inside map, filter, all, some and none`, node, text);
      return {var: ''};
    }
    return {var: path};
  }
  case 'Call': {
    const {name, positional, named} = callParts(node, text);
    if (has(effectCalls, name))
      throw new FormulaError(`${name} is an effect call and cannot be used inside an expression`, node, text);
    if (has(sourceCalls, name))
      throw new FormulaError(`${name} is a source call and cannot be used inside an expression`, node, text);
    if (name === 'obj') {
      if (positional.length)
        throw new FormulaError('obj takes named arguments only', positional[0], text);
      const value: Record<string, any> = {};
      for (const [key, arg] of Object.entries(named)) {
        const compiled = expression(arg, text, inElement);
        if (!isConstant(compiled))
          throw new FormulaError(`obj values must be constants (${key})`, arg, text);
        value[key] = unliteral(compiled);
      }
      return {literal: value};
    }
    const op = has(renamed, name) ? renamed[name] : name;
    // the second argument of an array op runs over each element, which reduce reads as current, not $it
    const scoped = (arg: IToken, idx: number) => {
      if (idx === 1 && (elementOps.has(name) || name === 'reduce'))
        return expression(arg, text, name !== 'reduce');
      return expression(arg, text, inElement);
    };
    if (!has(renamed, name) && !formulaOps.has(name))
      throw new FormulaError(`Unknown operation ${name}`, node, text);
    const firstNamed = Object.values(named)[0];
    if (firstNamed)
      throw new FormulaError(`${name} takes positional arguments only`, firstNamed, text);
    let args: any[];
    if (name === 'missing')
      args = positional.map((arg) => pathArg(arg, text, name));
    else if (name === 'missing_some') {
      const list = positional[1];
      if (positional.length !== 2 || list.type !== 'Array')
        throw new FormulaError('missing_some takes a count and a list of aliases', node, text);
      args = [expression(positional[0], text, inElement), list.children.map((arg) => pathArg(arg, text, name))];
    } else if (name === 'var') {
      if (positional.length < 1 || positional.length > 2)
        throw new FormulaError('var takes an alias and an optional default', node, text);
      const fallback = positional.slice(1).map((arg) => expression(arg, text, inElement));
      args = [pathArg(positional[0], text, name), ...fallback];
    } else
      args = positional.map(scoped);
    return {[op]: args};
  }
  default:
    throw new FormulaError(`Unexpected ${node.type}`, node, text);
  }
}

function targets(node: IToken, text: string): string | string[] {
  if (node.type === 'Ref' && !node.text.includes('.'))
    return node.text;
  if (node.type === 'Array' && node.children.every((child) => child.type === 'Ref' && !child.text.includes('.')))
    return node.children.map((child) => child.text);
  throw new FormulaError('Targets are an output alias or a list of them', node, text);
}

function alias(node: IToken, text: string, what: string): string {
  if (node.type !== 'Ref' || node.text.includes('.'))
    throw new FormulaError(`${what} is an alias`, node, text);
  return node.text;
}

function stringLiteral(node: IToken, text: string, what: string): string {
  if (node.type !== 'String')
    throw new FormulaError(`${what} must be a string`, node, text);
  return unquote(node, text);
}

function constant(node: IToken, text: string, what: string): any {
  const value = expression(node, text);
  if (!isConstant(value))
    throw new FormulaError(`${what} must be a constant`, node, text);
  return unliteral(value);
}

function rootCall(text: string, namespace: Record<string, unknown>, kind: string) {
  const root = parse(text);
  if (root.type !== 'Call' || !has(namespace, root.children[0].text)) {
    const name = root.type === 'Call' ? root.children[0].text : '';
    const got = root.type !== 'Call' ? (root.type === 'Ref' ? 'a reference' : 'a literal') :
      has(effectCalls, name) ? `an effect call ${name}` : has(sourceCalls, name) ? `a source call ${name}` :
        isOp(name) ? `an expression op ${name}` : `an unknown call ${name}`;
    throw new FormulaError(`Expected ${kind} call at the root, got ${got}`, root, text);
  }
  return {root, ...callParts(root, text)};
}

function arity(name: string, positional: IToken[], params: readonly string[], root: IToken, text: string) {
  if (positional.length > params.length)
    throw new FormulaError(`${name} takes ${params.length} positional argument(s) (${params.join(', ')})`, root, text);
  if (positional.length < params.length) {
    throw new FormulaError(
      `${name} needs ${params[positional.length]} (positional arguments: ${params.join(', ')})`, root, text);
  }
}

const isOp = (name: string) => has(renamed, name) || formulaOps.has(name) || name === 'obj';

const stripMarker = (text: string) => text.startsWith('=') ? text.slice(1) : text;

const condition = (prefix: string, when: RuleLogic | undefined): RuleLogic | undefined =>
  typeof when === 'string' ? at(prefix, () => compileExpression(when) as RuleLogic) : when;

/** A JSON Logic expression from formula text; a leading `=` is ignored. */
export function compileExpression(text: string): RuleExpr {
  const body = stripMarker(text);
  return expression(parse(body), body);
}

/** An effect object from an effect call such as `set(t, m, restriction: "restricted")`. */
export function compileEffect(text: string): RuleEffect {
  const body = stripMarker(text);
  const {root, name, positional, named} = rootCall(body, effectCalls, 'an effect');
  const {positional: params, options}: {positional: readonly string[], options: readonly string[]} =
    effectCalls[name as EffectName];
  arity(name, positional, params, root, body);
  const effect: Record<string, any> = {effect: name};
  params.forEach((param, idx) => {
    const arg = positional[idx];
    if (param === 'targets')
      effect.targets = targets(arg, body);
    else if (param === 'source')
      effect.source = alias(arg, body, 'source');
    else
      effect[param] = expression(arg, body);
  });
  if (name === 'meta')
    effect.meta = {};
  for (const [key, arg] of Object.entries(named)) {
    if (key === 'when')
      effect.when = expression(arg, body);
    else if (name === 'meta')
      effect.meta[key] = expression(arg, body);
    else if (!options.includes(key))
      throw new FormulaError(`${name} has no option ${key} (expected ${[...options, 'when'].join(', ')})`, arg, body);
    else if (key === 'targets')
      effect.targets = targets(arg, body);
    else if (key === 'restriction') {
      const restriction = stringLiteral(arg, body, 'restriction');
      if (!has(restrictions, restriction))
        throw new FormulaError(`Unknown restriction ${restriction}`, arg, body);
      effect.restriction = restriction;
    } else if (key === 'ignoreCase') {
      const value = constant(arg, body, 'ignoreCase');
      if (typeof value !== 'boolean')
        throw new FormulaError('ignoreCase must be true or false', arg, body);
      effect.ignoreCase = value;
    }
  }
  return effect as RuleEffect;
}

/** A source object from a source call such as `func("Pkg:F", key: key)`. */
export function compileSource(text: string): RuleSource {
  const body = stripMarker(text);
  const {root, name, positional, named} = rootCall(body, sourceCalls, 'a source');
  arity(name, positional, sourceCalls[name as keyof typeof sourceCalls].positional, root, body);
  const entries = Object.entries(named);
  const args = () => entries.length ?
    {args: Object.fromEntries(entries.map(([key, arg]) => [key, expression(arg, body)]))} : {};
  const noOptions = () => {
    if (entries.length)
      throw new FormulaError(`${name} has no option ${entries[0][0]}`, entries[0][1], body);
  };
  switch (name) {
  case 'file':
    noOptions();
    return {file: stringLiteral(positional[0], body, 'path')};
  case 'table': {
    const csv = stringLiteral(positional[0], body, 'csv');
    if (!entries.length)
      return {table: csv};
    const options = Object.fromEntries(entries.map(([key, arg]) => [key, constant(arg, body, key)]));
    return {table: {csv, options}};
  }
  case 'func':
    return {func: {name: stringLiteral(positional[0], body, 'name'), ...args()}};
  case 'query': {
    const connection = stringLiteral(positional[0], body, 'connection');
    return {query: {connection, sql: stringLiteral(positional[1], body, 'sql'), ...args()}};
  }
  case 'validators': {
    const input = alias(positional[0], body, 'input');
    const extra = entries.find(([key]) => key !== 'names');
    if (extra)
      throw new FormulaError(`validators has no option ${extra[0]} (expected names)`, extra[1], body);
    if (!named.names)
      return {validators: {input}};
    const names = constant(named.names, body, 'names');
    if (!Array.isArray(names) || names.some((item) => typeof item !== 'string'))
      throw new FormulaError('names must be a list of function names', named.names, body);
    return {validators: {input, names}};
  }
  default:
    noOptions();
    return {choices: {input: alias(positional[0], body, 'input')}};
  }
}

/** A value field: a string starting with `=` is a formula, any other value is kept as it is. */
export function compileValue(value: RuleExpr): RuleExpr {
  return typeof value === 'string' && value.startsWith('=') ? compileExpression(value) : value;
}

function at<T>(prefix: string, compile: () => T): T {
  try {
    return compile();
  } catch (e) {
    throw new Error(`${prefix}: ${e instanceof Error ? e.message : String(e)}`);
  }
}

const valueField = (prefix: string, value: RuleExpr): RuleExpr => at(prefix, () => compileValue(value));

const effectValueFields = ['items', 'message', 'value', 'values'] as const;

// malformed shapes are left as they are for expandRule, which reports them
const isMap = (value: unknown): value is Record<string, any> =>
  value != null && typeof value === 'object' && !Array.isArray(value);

function objectEffect(prefix: string, effect: RuleEffect): RuleEffect {
  if (!isMap(effect))
    return effect;
  const out: Record<string, any> = {...effect};
  out.when = condition(`${prefix}.when`, effect.when);
  for (const field of effectValueFields) {
    if (field in out)
      out[field] = valueField(`${prefix}.${field}`, out[field]);
  }
  if (effect.effect === 'meta' && isMap(effect.meta)) {
    out.meta = Object.fromEntries(Object.entries(effect.meta).map(([key, expr]) =>
      [key, valueField(`${prefix}.meta.${key}`, expr)]));
  }
  return out as RuleEffect;
}

function objectSource(prefix: string, source: RuleSource): RuleSource {
  const compileArgs = (args: Record<string, RuleExpr>) => Object.fromEntries(
    Object.entries(args).map(([key, expr]) => [key, valueField(`${prefix}.args.${key}`, expr)]));
  if (!isMap(source))
    return source;
  if ('func' in source && isMap(source.func?.args))
    return {func: {...source.func, args: compileArgs(source.func.args)}};
  if ('query' in source && isMap(source.query?.args))
    return {query: {...source.query, args: compileArgs(source.query.args)}};
  return source;
}

/** The rule with every formula string replaced by the object it stands for. */
export function compileRuleFormulas(rule: PipelineRuleConfiguration<any>): PipelineRuleConfiguration<any> {
  const prefix = `Rule ${rule.id}`;
  const when = condition(`${prefix}: when`, rule.when);
  const effects = !Array.isArray(rule.effects) ? rule.effects : rule.effects.map((effect, idx) =>
    typeof effect === 'string' ?
      at(`${prefix}: effects[${idx}]`, () => compileEffect(effect)) :
      objectEffect(`${prefix}: effects[${idx}]`, effect));
  const sources = !isMap(rule.sources) ? rule.sources : Object.fromEntries(
    Object.entries(rule.sources).map(([name, source]) => [name, typeof source === 'string' ?
      at(`${prefix}: sources.${name}`, () => compileSource(source)) :
      objectSource(`${prefix}: sources.${name}`, source)]));
  return {...rule, when, effects, ...(sources ? {sources} : {})};
}

/** The check with its `when` and `message` formulas compiled. */
export function compileCheckFormulas(check: PipelineCheckConfiguration<any>): PipelineCheckConfiguration<any> {
  const prefix = `Check ${check.id}`;
  const when = condition(`${prefix}: when`, check.when);
  const message = check.message === undefined ? undefined : valueField(`${prefix}: message`, check.message);
  return {...check, when, message};
}
