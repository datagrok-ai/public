//import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
//import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from './package';

export function getDivWithLatexContent(latexText: string): HTMLDivElement {
  const generator = new HtmlGenerator({hyphenate: false});
  const htmlGenerator = parse(latexText, {generator: generator});
  const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`);

  const div = ui.div([]);
  div.style.width = '100%';
  div.style.height = '100%';
  const iframe = ui.iframe({width: '100%', height: '100%'});
  div.append(iframe);

  iframe.srcdoc = doc.documentElement.outerHTML;

  return div;
}
