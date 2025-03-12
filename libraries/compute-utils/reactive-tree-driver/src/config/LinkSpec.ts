import {Grammars, Parser, IToken} from 'ebnf';
import {ItemId} from '../data/common-types';
import {indexFromEnd} from '../utils';
import {IOType} from './config-processing-utils';

const linkSpecGrammar = `
Link ::= Name FlagList? (':' Segment)? ('/' Segment)*
FlagList ::= WS* '(' WS* Flag WS* (',' WS* Flag WS*)* ')' WS* {fragment=true}
Flag ::= "call" | "optional"
Segment ::=  WS* (Selector | TargetIds | TagSpec) WS* {fragment=true}
TagSpec ::=  WS* '#' WS* TagSelectorType '(' TagSelectorArg ')' WS*
TagSelectorArg ::= ((RefArg ',' TagArg) | TagArg) {fragment=true}
TagSelectorType ::= "after*" | "after" | "before*" | "before" | "same" | "first" | "last" | "all" | "expand"
Selector ::= SelectorType '(' SelectorArgs ')'
SelectorArgs ::= ((RefArg ',' TargetIds) | TargetIds) (',' StopIds)? {fragment=true}
SelectorType ::= "after+" | "after*" | "after" | "before+" | "before*" | "before" | "same" | "first" | "last" | "all" | "expand"
TargetIds ::= ListArg
StopIds ::= ListArg
ListArg ::= WS* Id WS* ('|' WS* Id WS*)* {fragment=true}
RefArg ::= WS* '@' Ref WS* {fragment=true}
TagArg ::= WS* Tag WS* {fragment=true}
Ref ::= IDENTIFIER
Id ::= IDENTIFIER
Name ::= IDENTIFIER
Tag ::= IDENTIFIER
IDENTIFIER ::= [a-zA-Z][a-zA-Z_0-9]*
WS ::= ' '
`;

const linkParser = new Parser(Grammars.Custom.getRules(linkSpecGrammar));

export type LinkRefSelectors = 'after+' | 'after*' | 'after' | 'before+' | 'before*' | 'before' | 'same';
export type TagRefSelectors = 'after*' | 'after' | 'before*' | 'before' | 'same';
export type LinkNonRefSelectors = 'first' | 'last' | 'all' | 'expand';
export type LinkFlags = 'call' | 'optional';
export type LinkSelectors = LinkRefSelectors | LinkNonRefSelectors;
export type TagSelectors = LinkNonRefSelectors | TagRefSelectors;

export type LinkSelectorSegment = {
  type: 'selector',
  selector: LinkSelectors,
  ids: ItemId[],
  stopIds: ItemId[],
  ref?: string | undefined,
}

export type LinkTagSegment = {
  type: 'tag',
  selector: TagSelectors,
  tag: string,
  ref?: string | undefined,
}

export type LinkIOParsed = {
  name: string;
  segments: (LinkSelectorSegment | LinkTagSegment)[];
  flags?: LinkFlags[] | undefined,
}

const nonRefSelectors = ['first', 'last', 'all', 'expand'];

export function isNonRefSelector(sel: LinkSelectors): sel is LinkNonRefSelectors {
  return nonRefSelectors.includes(sel);
}

export function isRefSelector(sel: LinkSelectors): sel is LinkRefSelectors {
  return !nonRefSelectors.includes(sel);
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

export function parseLinkIO(io: string, ioType: IOType): LinkIOParsed {
  const ast = linkParser.getAST(io);
  checkAST(io, ast);
  const name = ast.children.find((cnode) => cnode.type === 'Name')!.text;
  const flags = ast.children.filter((cnode) => cnode.type === 'Flag').map((cnode) => cnode.text as LinkFlags);
  const isBase = ioType === 'base';
  const isAction = ioType === 'actions';
  const segments = ast.children.map((node) => {
    if (node.type === 'Selector') {
      const selector = node.children[0].text as LinkSelectors;
      const ref = node.children.find((cnode) => cnode.type === 'Ref')?.text;
      const targetIdsNode = node.children.find((cnode) => cnode.type === 'TargetIds');
      const stopIdsNode = node.children.find((cnode) => cnode.type === 'StopIds');
      checkSelector(io, selector, isBase, ref);
      const ids = targetIdsNode!.children.map((cnode) => cnode.text);
      const stopIds = stopIdsNode ? stopIdsNode.children.map((cnode) => cnode.text) : [];
      return {type: 'selector' as const, selector, ids, stopIds, ref};
    } else if (node.type === 'TargetIds') {
      const selector = 'first' as const;
      const ids = node.children.map((cnode) => cnode.text);
      return {type: 'selector' as const, ids, selector, stopIds: []};
    } else if(node.type === 'TagSpec') {
      const selector = node.children[0].text as TagSelectors;
      const ref = node.children.find((cnode) => cnode.type === 'Ref')?.text;
      const tag = node.children.find((cnode) => cnode.type === 'Tag')!.text;
      checkSelector(io, selector, isBase, ref);
      return {type: 'tag' as const, selector, tag, ref };
    } else if (node.type === 'Name' || node.type === 'Flag')
      return;
    throw new Error(`Link ${io}, unknown AST node type ${node.type}`);
  }).filter((x) => !!x);
  const lastSegment = indexFromEnd(segments);
  if (lastSegment && !isBase && !isAction && (lastSegment.selector !== 'first'))
    throw new Error(`Link io ${io} is not ending with input/output selector`);
  return {name, segments, flags};
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
