import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import JSZip from 'jszip';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';

const previews: [string, (file: DG.FileInfo) => Promise<DG.View>][] = [
  ['en-sdg.docx', PackageFunctions.previewDocx],
  ['en-sdg.doc', PackageFunctions.previewDoc],
  ['en-sdg.odt', PackageFunctions.previewOdt],
  ['es-sdg.odt', PackageFunctions.previewOdt],
  ['en-sdg.rtf', PackageFunctions.previewRtf],
];

const demoFile = async (name: string) =>
  (await grok.dapi.files.list('System:DemoFiles/texts/doc')).find((f) => f.name === name)!;

/** The rendered document: the preview itself, or the shadow root a DOCX renders into. */
const content = (view: DG.View): ParentNode => {
  const host = view.root.querySelector('.fe-document-preview')!;
  return host.querySelector('.fe-docx')?.shadowRoot ?? host;
};

const documentText = (view: DG.View) => Array.from(content(view).childNodes)
  .filter((n) => n.nodeName !== 'STYLE').map((n) => n.textContent).join('');

const links = (view: DG.View) => Array.from(content(view).querySelectorAll('a')).map((a) => a.getAttribute('href')).join();

/** Rewrites files inside a ZIP-based document (DOCX, ODT). */
async function patchZip(name: string, patches: Record<string, (xml: string) => string>): Promise<DG.FileInfo> {
  const zip = await JSZip.loadAsync(await (await demoFile(name)).readAsBytes());
  for (const [path, patch] of Object.entries(patches))
    zip.file(path, patch(await zip.file(path)!.async('string')));
  return DG.FileInfo.fromBytes(name, await zip.generateAsync({type: 'uint8array'}));
}

category('Documents', () => {
  for (const [name, preview] of previews) {
    test(name, async () => {
      const text = documentText(await preview(await demoFile(name)));
      expect(text.includes('cannot be displayed'), false, text.slice(0, 200));
      expect(text.length > 2000, true, `${name}: only ${text.length} characters`);
      expect(text.includes(name.startsWith('es') ? 'Objetivos de desarrollo sostenible' : 'Transforming our world'), true);
    }, {timeout: 60000});
  }

  test('Not a document', async () => {
    const text = documentText(await PackageFunctions.previewDoc(DG.FileInfo.fromString('x.doc', 'plain text')));
    expect(text.includes('cannot be displayed'), true, text);
  });

  test('DOCX symbol markup is not parsed as HTML', async () => {
    const file = await patchZip('es-sdg.docx', {'word/document.xml': (xml) => xml.replace('<w:r>',
      '<w:r><w:sym w:font="Symbol" w:char="41;&lt;img src=x onerror=void(0)&gt;"/>')});
    const view = await PackageFunctions.previewDocx(file);
    expect(content(view).querySelector('img[src="x"]'), null);
    expect(documentText(view).includes('Objetivos de desarrollo sostenible'), true);
  });

  test('DOCX styles stay inside the preview', async () => {
    const before = document.querySelectorAll('style').length;
    const view = await PackageFunctions.previewDocx(await demoFile('es-sdg.docx'));
    document.body.append(view.root);
    try {
      expect(document.querySelectorAll('style').length, before);
      expect(content(view).querySelectorAll('style').length > 0, true);
    }
    finally {
      view.root.remove();
    }
  });

  test('DOCX script links are dropped', async () => {
    const hyperlink = (target: string, id: string) =>
      `<Relationship Id="${id}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" ` +
      `Target="${target}" TargetMode="External"/>`;
    const file = await patchZip('es-sdg.docx', {
      'word/_rels/document.xml.rels': (xml) => xml.replace('</Relationships>',
        hyperlink('javascript:void(0)', 'feBad') + hyperlink('https://datagrok.ai', 'feGood') + '</Relationships>'),
      'word/document.xml': (xml) => xml.replace('<w:r>',
        '<w:hyperlink r:id="feBad"><w:r><w:t>bad</w:t></w:r></w:hyperlink>' +
        '<w:hyperlink r:id="feGood"><w:r><w:t>good</w:t></w:r></w:hyperlink><w:r>'),
    });
    expect(links(await PackageFunctions.previewDocx(file)), ',https://datagrok.ai');
  });

  test('ODT script links are dropped', async () => {
    const file = await patchZip('es-sdg.odt', {'content.xml': (xml) => xml.replace('<text:p ',
      '<text:p><text:a xlink:href="javascript:void(0)">bad</text:a><text:a xlink:href="https://datagrok.ai">good</text:a>' +
      '</text:p><text:p ')});
    expect(links(await PackageFunctions.previewOdt(file)), ',https://datagrok.ai');
  });

  test('RTF script links are dropped', async () => {
    const field = (url: string, text: string) => `{\\field{\\*\\fldinst HYPERLINK "${url}"}{\\fldrslt ${text}}}`;
    const rtf = `{\\rtf1\\ansi ${field('javascript:void(0)', 'bad')} ${field('https://datagrok.ai', 'good')}}`;
    expect(links(await PackageFunctions.previewRtf(DG.FileInfo.fromString('x.rtf', rtf))), ',https://datagrok.ai');
  });
});
