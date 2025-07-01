import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parse, HtmlGenerator} from 'latex.js';
import '../css/latex/css/katex.css';
import '../css/latex/css/article.css';


export function testLatex(): void {
  //const {parse, HtmlGenerator} = require('latex.js');

  //const latex = 'Hi, this is a line of text.';
  //

  const dlg = ui.dialog('Open a file');
  const fileInp = document.createElement('input');
  fileInp.type = 'file';
  fileInp.onchange = () => {
    //@ts-ignore
    const [file] = document.querySelector('input[type=file]').files;
    const reader = new FileReader();
    reader.addEventListener('load', () => {
      const latex = reader.result as string;

      const generator = new HtmlGenerator({hyphenate: false});

      const doc = parse(latex, {generator: generator}).htmlDocument();

      console.log(doc.documentElement.outerHTML);

      //const generator = new HtmlGenerator({hyphenate: false});

      //const doc = parse(latex, {generator: generator}).htmlDocument();

      console.log(doc);

      console.log(doc.documentElement.outerHTML);

      const div = ui.div([]);

      const v = grok.shell.addView(DG.View.create());

      v.append(div);

      div.innerHTML = doc.documentElement.outerHTML;

      dlg.close();
    }, false);

    if (file)
      reader.readAsText(file);
  };

  dlg.add(fileInp);
  fileInp.click();

//   ui.dialog({title: 'LaTex'})
//     .add(div)
//     .show();
}

export async function texFilePreview(file: DG.FileInfo): Promise<DG.View> {
  const latexText = await file.readAsString();
  const generator = new HtmlGenerator({hyphenate: false});
  const doc = parse(latexText, {generator: generator}).htmlDocument();
  const div = ui.div([]);
  const v = DG.View.create();
  v.append(div);
  v.name = file.fileName;
  div.innerHTML = doc.documentElement.outerHTML;

  return v;
} // texFilePreview
