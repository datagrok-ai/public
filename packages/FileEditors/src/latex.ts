import * as ui from 'datagrok-api/ui';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from './package';

export function getElementWithLatexContent(latexText: string): HTMLIFrameElement {
  const generator = new HtmlGenerator({hyphenate: false});
  const htmlGenerator = parse(latexText, {generator: generator});
  const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`) as Document;

  // const div = ui.div([]);
  // div.style.width = '100%';
  // div.style.height = '100%';
  const iframe = ui.iframe({width: '100%', height: '100%'});
  //div.append(iframe);

  const rightMargin = doc.body.children[1];

  if ((rightMargin !== null) && (rightMargin !== null))
    (rightMargin as HTMLElement).style.display = 'none';

  iframe.srcdoc = doc.documentElement.outerHTML;
  iframe.style.border = 'transparent';

  return iframe;
}
