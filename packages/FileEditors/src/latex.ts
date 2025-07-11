import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
//import '../css/latex/css/katex.css';
//import '../css/latex/css/article.css';

// import katex from './katex/katex.min.js';
// import './katex/katex.min.css';

export async function texFilePreview(file: DG.FileInfo): Promise<DG.View> {
  try {
    const latexText = await file.readAsString();
    const generator = new HtmlGenerator({hyphenate: false});
    const htmlGenerator = parse(latexText, {generator: generator});

    //console.log(htmlGenerator);

    const doc = htmlGenerator.htmlDocument();

    //console.log(doc);

    const div = ui.div([]);
    const v = DG.View.create();
    v.append(div);

    div.style.width = '100%';
    div.style.height = '100%';

    v.name = file.fileName;

    //console.log(doc.documentElement);

    //console.log(doc.documentElement.outerHTML);

    //div.innerHTML = doc.documentElement.outerHTML;
    //div.append(doc.documentElement);

    //const domFragment = htmlGenerator.domFragment();
    //console.log(domFragment);

    //div.append(domFragment);

    // Create and inject into iframe
    const html = doc.documentElement.outerHTML;
    console.log((html as string).indexOf('https://dev.datagrok.ai/css/katex.css'));
    console.log((html as string).indexOf('https://dev.datagrok.ai/css/article.css'));
    const iframe = ui.iframe({src: html, width: '100%', height: '100%'});
    div.append(iframe);

    // console.log('contentDocument:', iframe.contentDocument);
    // console.log('contentWindow:', iframe.contentWindow);

    // console.log('iframe:', iframe);

    iframe.srcdoc = html;

    // const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
    // iframeDoc.open();
    // iframeDoc.write(doc.documentElement.outerHTML);
    // iframeDoc.close();

    return v;
  } catch (err) {
    const v = DG.View.create();
    v.append(ui.div([ui.h2('The file is corrupted and cannot be opened!')]));
    v.name = file.fileName;

    if (err instanceof Error)
      grok.shell.error(err.message);

    return v;
  }
  // try {
  //   const latexText = await file.readAsString();
  //   const div = ui.div([]);

  //   katex.render(latexText, div, {
  //     throwOnError: false,
  //   });

  //   const v = DG.View.create();
  //   v.append(div);
  //   v.name = file.fileName;

  //   return v;
  // } catch (err) {
  //   const v = DG.View.create();
  //   v.append(ui.div([ui.h2('The file is corrupted and cannot be opened!')]));
  //   v.name = file.fileName;

  //   if (err instanceof Error)
  //     grok.shell.error(err.message);

  //   return v;
  // }
} // texFilePreview
