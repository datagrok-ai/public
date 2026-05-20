import {Grammars, Parser, IToken} from 'ebnf';
import {ItemId, NqName} from '../data/common-types';
import {indexFromEnd} from '../utils';
import {IOType} from './config-processing-utils';

/* eslint-disable max-len */
const linkSpecGrammar = `
Link ::= Name FlagList? (':' Tail?)?
Tail ::= '.' '/' Segments? | '.' '/'? | Segments {fragment=true}
Segments ::= (Segment '/')* (IOSegment | Segment) {fragment=true}
FlagList ::= WS* '(' WS* Flag WS* (',' WS* Flag WS*)* ')' WS* {fragment=true}
Flag ::= "call" | "optional" | "template"
Segment ::= WS* (Selector | TargetIds | TagSpec) WS* {fragment=true}
IOSegment ::= WS* IOSelector WS* {fragment=true}
TagSpec ::= WS* '#' WS* TagSelectorType '(' TagSelectorArgs ')' WS*
TagSelectorArgs ::= (RefArg (',' TagIds)?) | TagIds {fragment=true}
TagSelectorType ::= "after*" | "after" | "before*" | "before" | "same" | "first" | "last" | "all" | "expand"
IOSelector ::= IOSelectorType '(' WS* IOSelectorNqName WS* (',' WS* TargetIds?)? WS* ')'
IOSelectorType ::= "inputs" | "outputs"
IOSelectorNqName ::= Id ':' Id
Selector ::= SelectorType '(' SelectorArgs? ')'
SelectorArgs ::= (RefArg (',' TargetIds (',' StopIds)?)?) | (TargetIds (',' StopIds)?) {fragment=true}
SelectorType ::= "after+" | "after*" | "after" | "before+" | "before*" | "before" | "same" | "first" | "last" | "all" | "expand"
TargetIds ::= ArgsOr
StopIds ::= ArgsOr
TagIds ::= ArgsAnd
ArgsOr ::= WS* Id WS* ('|' WS* Id WS*)* {fragment=true}
ArgsAnd ::= WS* Id WS* ('&' WS* Id WS*)* {fragment=true}
RefArg ::= WS* '@' Ref WS* {fragment=true}
Ref ::= IDENTIFIER
Id ::= IDENTIFIER
Name ::= IDENTIFIER
IDENTIFIER ::= [_a-zA-Z][a-zA-Z_0-9]*
WS ::= ' '
`;
/* eslint-enable max-len */

const linkParser = new Parser(Grammars.Custom.getRules(linkSpecGrammar));

export type LinkRefSelectors = 'after+' | 'after*' | 'after' | 'before+' | 'before*' | 'before' | 'same';
export type TagRefSelectors = 'after*' | 'after' | 'before*' | 'before' | 'same';
export type LinkNonRefSelectors = 'first' | 'last' | 'all' | 'expand';
export type IOExpandSelectors = 'inputs' | 'outputs';
export type LinkFlags = 'call' | 'optional' | 'template';
export type LinkSelectors = LinkRefSelectors | LinkNonRefSelectors;
export type TagSelectors = LinkNonRefSelectors | TagRefSelectors;

export type LinkSelectorSegment = {
  type: 'selector',
  selector: LinkSelectors,
  ids: ItemId[],
  stopIds: ItemId[],
  ref?: string | undefined,
  // Set on the last segment of a `(template)` query that used `inputs(...)`
  // or `outputs(...)`. Consumed at config-processing time, where the
  // deferred entry is replaced with N concrete bare-list LinkIOParsed
  // entries. The runtime never sees this.
  ioExpand?: IOExpandSelectors,
  nqName?: NqName,
  excludeIds?: ItemId[],
}

export type LinkTagSegment = {
  type: 'tag',
  selector: TagSelectors,
  tags: string[],
  ref?: string | undefined,
}

export type LinkIOParsed = {
  name: string;
  segments: (LinkSelectorSegment | LinkTagSegment)[];
  flags?: LinkFlags[] | undefined,
  // Canonical slot id for the controller API (propagateTemplatePair,
  // getInputTemplates/getOutputTemplates). For prefixed operators
  // (e.g. `in_(template)`) it is the prefix string `'in_'`. For the
  // anonymous form `_(template)` it is the 0-based numeric index of the
  // anonymous operator on its side — using a number prevents any collision
  // with user-typed prefix strings (grammar identifiers cannot be pure digits).
  templateName?: string | number;
}

const nonRefSelectors = ['first', 'last', 'all', 'expand'];

export function isNonRefSelector(sel: LinkSelectors): sel is LinkNonRefSelectors {
  return nonRefSelectors.includes(sel);
}

export function isRefSelector(sel: LinkSelectors): sel is LinkRefSelectors {
  return !nonRefSelectors.includes(sel);
}

export function isDeferredIOSegment(seg: LinkSelectorSegment | LinkTagSegment): boolean {
  return seg.type === 'selector' && !!seg.ioExpand;
}

export function refSelectorAdjacent(sel: LinkRefSelectors) {
  return sel.endsWith('+');
}

export function refSelectorAll(sel: LinkRefSelectors) {
  return sel.endsWith('*');
}

export function refSelectorFindOne(sel: LinkRefSelectors) {
  return sel === 'after' || sel === 'before';
}

export type SelectorDirection = 'before' | 'after' | 'same';

export function refSelectorDirection(sel: LinkRefSelectors): SelectorDirection {
  if (sel === 'same')
    return 'same' as const;
  return sel.startsWith('after') ? 'after' : 'before';
}

export function parseLinkIO(io: string, ioType: IOType): LinkIOParsed[] {
  const ast = linkParser.getAST(io);
  checkAST(io, ast);
  const name = ast.children.find((cnode) => cnode.type === 'Name')!.text;
  const flags = ast.children.filter((cnode) => cnode.type === 'Flag').map((cnode) => cnode.text as LinkFlags);
  const isBase = ioType === 'base';
  const isAction = ioType === 'actions';
  const isNot = ioType === 'not';
  const isTemplate = flags.includes('template');
  if (flags.includes('call')) {
    if (isTemplate)
      throw new Error(`Link io ${io}: (call) and (template) flags cannot be combined`);
    if (isBase)
      throw new Error(`Link io ${io}: (call) flag is not allowed in base queries`);
  }

  type Segment = LinkSelectorSegment | LinkTagSegment;
  const segments = ast.children.map((node): Segment | undefined => {
    if (node.type === 'Selector') {
      const selector = node.children[0].text as LinkSelectors;
      const ref = node.children.find((cnode) => cnode.type === 'Ref')?.text;
      const targetIdsNode = node.children.find((cnode) => cnode.type === 'TargetIds');
      const stopIdsNode = node.children.find((cnode) => cnode.type === 'StopIds');
      checkSelector(io, selector, isBase, ref);
      const ids = targetIdsNode ? targetIdsNode.children.map((cnode) => cnode.text) : [];
      const stopIds = stopIdsNode ? stopIdsNode.children.map((cnode) => cnode.text) : [];
      if (ids.length === 0 && selector !== 'same')
        throw new Error(`Link io ${io} is using ${selector} without an id list`);
      return {type: 'selector' as const, selector, ids, stopIds, ref};
    } else if (node.type === 'IOSelector') {
      const ioExpand = node.children[0].text as IOExpandSelectors;
      const nqName = node.children.find((c) => c.type === 'IOSelectorNqName')!.text.trim();
      const targetIdsNode = node.children.find((c) => c.type === 'TargetIds');
      const excludeIds = targetIdsNode ? targetIdsNode.children.map((c) => c.text) : [];
      if (isBase || isAction || isNot)
        throw new Error(`Link io ${io}: '${ioExpand}' selector is not allowed in ${ioType} queries`);
      return {type: 'selector' as const, selector: 'first', ids: [], stopIds: [], ioExpand, nqName, excludeIds};
    } else if (node.type === 'TargetIds') {
      const ids = node.children.map((cnode) => cnode.text);
      return {type: 'selector' as const, ids, selector: 'first', stopIds: []};
    } else if (node.type === 'TagSpec') {
      const selector = node.children[0].text as TagSelectors;
      const ref = node.children.find((cnode) => cnode.type === 'Ref')?.text;
      const tagsNode = node.children.find((cnode) => cnode.type === 'TagIds');
      const tags = tagsNode ? tagsNode.children.map((cnode) => cnode.text) : [];
      checkSelector(io, selector, isBase, ref);
      if (tags.length === 0 && selector !== 'same')
        throw new Error(`Link io ${io} is using ${selector} without a tag list`);
      return {type: 'tag' as const, selector, tags, ref};
    } else if (node.type === 'Name' || node.type === 'Flag')
      return;
    throw new Error(`Link ${io}, unknown AST node type ${node.type}`);
  }).filter((x): x is Segment => !!x);

  const lastSeg = indexFromEnd(segments);
  const isDeferred = !!(lastSeg && lastSeg.type === 'selector' && lastSeg.ioExpand);

  // Grammar guarantees IO selectors only appear as the last segment. We still
  // need to enforce that they require the (template) flag.
  if (isDeferred && !isTemplate)
    throw new Error(`Link io ${io}: '${(lastSeg as LinkSelectorSegment).ioExpand}' selector requires the (template) flag`);

  if (lastSeg && !isBase && !isAction && lastSeg.type === 'selector' &&
      lastSeg.selector !== 'first' && !isDeferred)
    throw new Error(`Link io ${io} is not ending with input/output selector`);

  // Deferred IO: emit a single LinkIOParsed; expansion happens at config-time.
  if (isDeferred)
    return [{name, segments, flags}];

  // Existing bare-list template expansion (unchanged).
  if (lastSeg && lastSeg.type === 'selector' && isTemplate) {
    return lastSeg.ids!.map((id) => {
      const nname = name === '_' ? id : name + id;
      const nlastSegment: LinkSelectorSegment = {type: 'selector', selector: 'first', ids: [id], stopIds: []};
      const nsegments = [...segments.slice(0, -1), nlastSegment];
      return {name: nname, segments: nsegments, flags};
    });
  }

  return [{name, segments, flags}];
}

function checkSelector(io: string, selector: LinkSelectors | TagSelectors, isBase: boolean, ref?: string) {
  if (ref && isNonRefSelector(selector))
    throw new Error(`Link io ${io} is using non-ref selector with ref`);
  if (isBase && !isNonRefSelector(selector))
    throw new Error(`Link io ${io} is using ref selector ${selector} in link base`);
  if (isBase && selector === 'all')
    throw new Error(`Link io ${io} is using all selector in link base`);
  if (!isBase && selector === 'expand')
    throw new Error(`Link io ${io} is using expand selector in non-base`);
}

function checkAST(str: string, ast?: IToken) {
  if (ast == null)
    throw new Error(`Failed to parse link spec: ${str}`);
  if (ast.errors?.length)
    throw new Error(`Failed to parse link spec: ${str}, errors: ${ast.errors.map((e) => e.message).join(',')}`);
}
