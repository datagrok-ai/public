/** Documents on the program page (F-27): the dialog that links a file as a document (a file of the
 * program's workspace folder, any file path, or an upload into that folder), and the viewer a
 * document opens in. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {ProgramRow} from '../generated/db';
import {folderFiles, LinkedDocument, linkDocumentTo, uploadDocument} from '../documents';
import {describeDocument, DocumentKindConfig, DocumentStatus} from '../model/documents';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

function field(label: string, control: HTMLElement, name: string): HTMLElement {
  return named(divH([span(label, 'pp-field-label'), control], 'pp-field'), name);
}

function select(options: {code: string, name: string}[], value: string): HTMLSelectElement {
  const el = document.createElement('select');
  el.className = 'pp-select';
  el.append(...options.map((o) => Object.assign(document.createElement('option'), {value: o.code, textContent: o.name})));
  el.value = value;
  return el;
}

function text(value = '', placeholder = ''): HTMLInputElement {
  const el = document.createElement('input');
  el.className = 'pp-text';
  el.value = value;
  el.placeholder = placeholder;
  return el;
}

export interface LinkDocumentOptions {
  program: ProgramRow;
  /** The document criterion the card belongs to; null to link a document no criterion asks for. */
  criterion: {id: string, code: string, name: string} | null;
  vocabulary: DocumentKindConfig[];
  /** Approvers may link a document approved; everyone else links drafts. */
  approver: boolean;
  onLinked?: (r: LinkedDocument) => void | Promise<void>;
}

export async function openLinkDocumentDialog(opts: LinkDocumentOptions): Promise<DG.Dialog> {
  installStyles();
  const {program, criterion, vocabulary} = opts;
  const files = await folderFiles(program).catch(() => []);
  const folderFile = select([{code: '', name: files.length ? 'Choose a file…' : 'No files in the folder yet'},
    ...files.map((f) => ({code: f.path, name: f.name}))], '');
  const path = text('', 'Or a file path: System:AppData/…');
  const upload = named(document.createElement('input'), 'document-upload-input');
  upload.type = 'file';
  upload.disabled = !program.space_id;
  upload.title = program.space_id ? 'Upload into this program\'s workspace folder' : 'No workspace folder: link a file by its path';
  const kind = select(vocabulary.map((k) => ({code: k.code, name: k.name})),
    vocabulary.find((k) => k.code === criterion?.code)?.code ?? vocabulary[vocabulary.length - 1].code);
  const title = text('', criterion ? `${criterion.name}, ${program.name}` : 'Title');
  const version = text('', '1.0');
  const statuses: DocumentStatus[] = opts.approver ? ['draft', 'approved'] : ['draft'];
  const status = select(statuses.map((s) => ({code: s, name: s})), 'draft');
  const errors = named(divV([], 'pp-warnings'), 'document-errors');

  const dlg = ui.dialog({title: `Link a document${criterion ? `: ${criterion.name}` : ''}`});
  dlg.root.dataset.u2Name = 'link-document-dialog';
  const ok = () => dlg.getButton('LINK');
  const chosen = () => upload.files?.[0]?.name ?? (folderFile.value || path.value.trim());
  const refresh = () => {
    const problems = chosen() ? [] : ['Choose a file of the folder, type a path, or upload a file'];
    errors.replaceChildren(...problems.map((p) => span(p, 'pp-muted')));
    ok()?.classList.toggle('disabled', problems.length > 0);
    // a title from the file until one is typed
    if (!title.value && chosen())
      title.placeholder = chosen().split(/[/\\:]/).pop()!.replace(/\.[^.]+$/, '');
    return problems;
  };
  folderFile.onchange = refresh;
  path.oninput = refresh;
  upload.onchange = refresh;

  dlg.add(named(divV([
    span(`${program.code} · ${program.name}`, 'pp-muted'),
    field('File in the folder', folderFile, 'document-file-input'),
    field('File path', path, 'document-path-input'),
    field('Upload', upload, 'document-upload'),
    field('Kind', kind, 'document-kind-input'),
    field('Title', title, 'document-title-input'),
    field('Version', version, 'document-version-input'),
    field('Status', status, 'document-status-input'),
    ...(opts.approver ? [] : [span('Approvers approve: approving meets the criterion', 'pp-muted')]),
    errors,
  ], 'pp-save-run'), 'link-document-form'));

  let uploaded: {file: File, path: string, folder: string | null} | null = null;
  dlg.addButton('LINK', async () => {
    if (refresh().length > 0 || ok()?.classList.contains('pp-busy'))
      return;
    ok()?.classList.add('pp-busy', 'disabled');
    try {
      let target = folderFile.value || path.value.trim();
      let folder: string | null = folderFile.value ? program.space_id ?? null : null;
      const file = upload.files?.[0];
      if (file) {
        // uploaded once: a retry after a failed link links the file already in the folder
        if (uploaded?.file !== file)
          uploaded = {file, ...await uploadDocument(program, file.name, new Uint8Array(await file.arrayBuffer()))};
        target = uploaded.path;
        folder = uploaded.folder;
      }
      const r = await linkDocumentTo(program, {path: target, folder, criterionId: criterion?.id ?? null, kind: kind.value,
        title: title.value.trim() || title.placeholder, version: version.value.trim() || undefined,
        status: status.value as DocumentStatus});
      dlg.close();
      grok.shell.info([`${program.code}: ${describeDocument({title: r.document.title, version: r.document.doc_version,
        status: r.document.status as DocumentStatus})} linked` + (r.criterion ? `; ${r.criterion.code} ${r.criterion.state.replace('_', ' ')}` : ''),
      ...r.notes].join('\n'));
      await opts.onLinked?.(r);
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'link-document-error'));
      ok()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 560});
  named(ok(), 'link-document-button');
  refresh();
  return dlg;
}

const TEXT_FILE = /\.(md|txt|json|ya?ml|xml|log)$/i;

/** Opens a document's file: a table for CSV, the text in a dialog for text files, else a download. */
export async function openDocumentFile(d: {title: string, path: string, doc_version?: string | null, status?: string}):
  Promise<void> {
  if (/\.(csv|tsv)$/i.test(d.path)) {
    grok.shell.addTableView(await grok.dapi.files.readCsv(d.path));
    return;
  }
  const name = d.path.split(/[/\\:]/).pop() ?? d.path;
  if (TEXT_FILE.test(d.path)) {
    const body = await grok.dapi.files.readAsText(d.path);
    const dlg = ui.dialog({title: describeDocument({title: d.title, version: d.doc_version, status: d.status as DocumentStatus})});
    dlg.root.dataset.u2Name = 'document-viewer';
    const pre = named(document.createElement('pre'), 'document-text');
    pre.className = 'pp-card-text';
    pre.textContent = body;
    dlg.add(divV([span(d.path, 'pp-muted'), pre]));
    dlg.show({width: 640});
    return;
  }
  DG.Utils.download(name, new Blob([(await grok.dapi.files.readAsBytes(d.path)).slice()]));
}
