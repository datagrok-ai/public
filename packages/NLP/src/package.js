/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: supportedExt
//input: string filename
//output: bool supported
export function supportedExt(filename) {
    let extensions = `.csv, .doc, .docx, .eml, .epub, .gif, .htm, .html,
    .jpeg, .jpg, .json, .log, .mp3, .msg, .odt, .ogg, .pdf, .png, .pptx, .ps, .psv,
    .rtf, .tff, .tif, .tiff, .tsv, .txt, .wav, .xls, .xlsx`.replace(/\s{4}/g, '');
    return extensions.split(', ').some(ext => filename.endsWith(ext));
}
