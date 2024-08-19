import { Grammars, IToken } from 'ebnf';
import {ItemId} from '../data/common-types';

const fromGrammar = `
From ::= (Name ':')? (Selector | Id) (SEPARATOR (Selector | Id))*
Selector ::= ( SelectorType '(' WS* Id WS* ')' )
SelectorType ::= "all" | "last" | "first"
Id ::= IDENTIFIER
Name ::= IDENTIFIER
IDENTIFIER ::= [a-zA-Z][a-zA-Z_0-9]*
SEPARATOR ::= '/'
WS ::= ' '
`;

const fromParser = new Grammars.W3C.Parser(fromGrammar);

const toGrammar = `
To ::= (Name ':')? (Selector | Id) (SEPARATOR (Selector | Id))*
Selector ::= ( SelectorType '(' WS* Id WS* (',' WS* '@' Ref WS*)? SelectorArgs* ')' )
SelectorArgs ::= (',' WS* Arg WS*)*
SelectorType ::= "nextUntil" | "next" | "adjacent" | "allUntil" | "all"
Id ::= IDENTIFIER
Name ::= IDENTIFIER
Arg ::= IDENTIFIER
Ref ::= IDENTIFIER
IDENTIFIER ::= [a-zA-Z][a-zA-Z_0-9]*
SEPARATOR ::= '/'
WS ::= ' '
`;

const toParser = new Grammars.W3C.Parser(toGrammar);

export type InputSelectors = 'last' | 'all' | 'first';

export type LinkInputSegment = {
  id: ItemId,
  selector: InputSelectors,
}

export type LinkInput = {
  name?: string;
  segments: LinkInputSegment[];
}

export type OutputSelectors = 'next' | 'nextUntil' | 'adjacent' | 'all' | 'allUntil';

export type LinkOutputSegment = {
  id: ItemId,
  selector: OutputSelectors,
  args?: string[],
  ref?: string | undefined,
}

export type LinkOutput = {
  name?: string;
  segments: LinkOutputSegment[];
}


export function parseLinkInput(input: string): LinkInput {
  const ast = fromParser.getAST(input);
  checkAST(input, ast);
  const name = ast.children.find(cnode => cnode.type === 'Name')?.text;
  const segments = ast.children.map((node) => {
    if (node.type === 'Selector') {
      const selector = node.children[0].text as InputSelectors;
      const id = node.children[1].text;
      return {id, selector};
    } else if (node.type === 'Id') {
      const id = node.text;
      const selector = 'last' as const;
      return {id, selector};
    } else if (node.type === 'Name') {
      return;
    }
    throw new Error(`Link input ${input}, unknown AST node type ${node.type}`);
  }).filter(x => !!x);
  return {name, segments};
}

export function parseLinkOutput(output: string): LinkOutput {
  const ast = toParser.getAST(output);
  checkAST(output, ast);
  const name = ast.children.find(cnode => cnode.type === 'Name')?.text;
  const segments = ast.children.map((node) => {
    if (node.type === 'Selector') {
      const selector = node.children[0].text as OutputSelectors;
      const id = node.children[1].text;
      const refNode = node.children.find((cnode) => cnode.type === 'Ref');
      const argsNode = node.children.find((cnode) => cnode.type === 'SelectorArgs');
      const ref = refNode ? refNode.text : undefined;
      const args = argsNode ? argsNode.children.map((cnode) => cnode.text) : [];
      return {id, selector, ref, args};
    } else if (node.type === 'Id') {
      const id = node.text;
      const selector = 'next' as const;
      return {id, selector};
    } else if (node.type === 'Name') {
      return;
    }
    throw new Error(`Link output ${output}, unknown AST node type ${node.type}`);
  }).filter(x => !!x);
  return {name, segments};
}

function checkAST(str: string, ast?: IToken) {
  if (ast == null)
    throw new Error(`Failed to parse link spec: ${str}`);
  if (ast.errors?.length)
    throw new  Error(`Failed to parse link spec: ${str}, errors: ${ast.errors.map(e => e.message).join(',')}`);
}
