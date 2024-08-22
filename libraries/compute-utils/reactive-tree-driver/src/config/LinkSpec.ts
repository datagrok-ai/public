import {Grammars, IToken} from 'ebnf';
import {ItemId} from '../data/common-types';

const linkSpecGrammar = `
To ::= Name ':' (Selector | Id) (SEPARATOR (Selector | Id))*
Selector ::= ( SelectorType '(' WS* Id WS* (',' WS* '@' Ref WS*)? SelectorArgs* ')' )
SelectorArgs ::= (',' WS* Arg WS*)*
SelectorType ::= "after+" | "after*" | "after" | "before+" | "before*" | "before" | "same" | "first" | "last" | "all" | "any"
Id ::= IDENTIFIER
Name ::= IDENTIFIER
Arg ::= IDENTIFIER
Ref ::= IDENTIFIER
IDENTIFIER ::= [a-zA-Z][a-zA-Z_0-9]*
SEPARATOR ::= '/'
WS ::= ' '
`;

const linkParser = new Grammars.W3C.Parser(linkSpecGrammar);

export type LinkRefSelectors = 'after+' | 'after*' | 'after' | 'before+' | 'before*' | 'before' | 'same';
export type LinkNonRefSelectors = 'first' | 'last' | 'all' | 'any';
export type LinkSelectors = LinkRefSelectors | LinkNonRefSelectors;

export type LinkSegment = {
  id: ItemId,
  selector: LinkSelectors,
  args?: string[],
  ref?: string | undefined,
}

export type LinkParsed = {
  name: string;
  segments: LinkSegment[];
}

const nonRefSelectors = ['first', 'last', 'all', 'any'];

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

export function parseLinkIO(io: string, isBase: boolean, isInput: boolean): LinkParsed {
  const ast = linkParser.getAST(io);
  checkAST(io, ast);
  const name = ast.children.find((cnode) => cnode.type === 'Name')!.text;
  const segments = ast.children.map((node) => {
    if (node.type === 'Selector') {
      const selector = node.children[0].text as LinkSelectors;
      const id = node.children[1].text;
      const refNode = node.children.find((cnode) => cnode.type === 'Ref');
      const argsNode = node.children.find((cnode) => cnode.type === 'SelectorArgs');
      const ref = refNode ? refNode.text : undefined;
      if (ref && isNonRefSelector(selector))
        throw new Error(`Link io ${io} is using non-ref selector with ref`);
      if (isBase && !isNonRefSelector(selector))
        throw new Error(`Link io ${io} is using ref selector ${selector} in link base`);
      if (isBase && selector === 'all')
        throw new Error(`Link io ${io} is using all selector in link base`);
      if (!isBase && selector === 'any')
        throw new Error(`Link io ${io} is using any selector in non-base`);
      const args = argsNode ? argsNode.children.map((cnode) => cnode.text) : [];
      return {id, selector, ref, args};
    } else if (node.type === 'Id') {
      const id = node.text;
      const selector = isInput ? ('last' as const) : ('first' as const);
      return {id, selector};
    } else if (node.type === 'Name')
      return;

    throw new Error(`Link ${io}, unknown AST node type ${node.type}`);
  }).filter((x) => !!x);
  return {name, segments};
}

function checkAST(str: string, ast?: IToken) {
  if (ast == null)
    throw new Error(`Failed to parse link spec: ${str}`);
  if (ast.errors?.length)
    throw new Error(`Failed to parse link spec: ${str}, errors: ${ast.errors.map((e) => e.message).join(',')}`);
}
