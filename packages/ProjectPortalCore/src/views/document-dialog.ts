/** Documents on the program page (F-27): the dialog that links a file as a document (a file of the
 * program's workspace folder, any file path, or an upload into that folder), and the viewer a
 * document opens in. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {divH, divV, signal, span} from '@datagrok-libraries/u2';
import {ProgramRow} from '../generated/db';
import {folderFiles, LinkedDocument, linkDocumentTo, uploadDocument} from '../documents';
import {describeDocument, DocumentKindConfig, DocumentStatus} from '../model/documents';
import {choice, choiceValue, ErrorRegion, named, openDialog, PortalDialog, text} from './u2-controls';

/** The raw file input keeps the old label-and-control row (u2 has no plain local-file input). */
const field = (label: string, control: HTMLElement, name: string): HTMLElement =>
  named(divH([span(label, 'pp-field-label'), control], 'pp-field'), name);

export interface LinkDocumentOptions {
  program: ProgramRow;
  /** The document criterion the card belongs to; null to link a document no criterion asks for. */
  criterion: {id: string, code: string, name: string} | null;
  vocabulary: DocumentKindConfig[];
  /** Approvers may link a document approved; everyone else links drafts. */
  approver: boolean;
  onLinked?: (r: LinkedDocument) => void | Promise<void>;
}

export async function openLinkDocumentDialog(opts: LinkDocumentOptions): Promise<PortalDialog> {
  const {program, criterion, vocabulary} = opts;
  const files = await folderFiles(program).catch(() => []);
  const upload = named(document.createElement('input'), 'document-upload-input');
  upload.type = 'file';
  upload.disabled = !program.space_id;
  upload.title = program.space_id ? 'Upload into this program\'s workspace folder' :
    'No workspace folder: link a file by its path';
  const statuses: DocumentStatus[] = opts.approver ? ['draft', 'approved'] : ['draft'];
  const errors = new ErrorRegion('document-errors', 'link-document-error');
  const valid = signal(false);
  let uploaded: {file: File, path: string, folder: string | null} | null = null;
  let linked: LinkedDocument;

  return openDialog({title: `Link a document${criterion ? `: ${criterion.name}` : ''}`, name: 'link-document-dialog',
    width: 560, errors}, () => {
    const fileOptions = files.map((f) => ({code: f.path, name: f.name}));
    const folderFile = choice('document-file-input', 'File in the folder', fileOptions, '',
      {blank: files.length ? 'Choose a file…' : 'No files in the folder yet', onChanged: () => refresh()});
    const path = text('document-path-input', 'File path', {placeholder: 'Or a file path: System:AppData/…'});
    const kind = choice('document-kind-input', 'Kind', vocabulary.map((k) => ({code: k.code, name: k.name})),
      vocabulary.find((k) => k.code === criterion?.code)?.code ?? vocabulary[vocabulary.length - 1].code);
    const title = text('document-title-input', 'Title',
      {placeholder: criterion ? `${criterion.name}, ${program.name}` : 'Title'});
    const version = text('document-version-input', 'Version', {placeholder: '1.0'});
    const status = choice('document-status-input', 'Status', statuses.map((x) => ({code: x, name: x})), 'draft');
    const chosen = () => upload.files?.[0]?.name ?? (choiceValue(folderFile) || path.value.value.trim());
    const refresh = () => {
      const problems = chosen() ? [] : ['Choose a file of the folder, type a path, or upload a file'];
      errors.root.replaceChildren(...problems.map((p) => span(p, 'pp-muted')));
      valid.value = problems.length === 0;
      // a title from the file until one is typed
      if (!title.value.value && chosen())
        title.placeholder = chosen().split(/[/\\:]/).pop()!.replace(/\.[^.]+$/, '');
    };
    path.effect(() => {
      void path.value.value;
      refresh();
    });
    upload.onchange = refresh;
    return {
      content: named(divV([
        span(`${program.code} · ${program.name}`, 'pp-muted'),
        folderFile.root, path.root,
        field('Upload', upload, 'document-upload'),
        kind.root, title.root, version.root, status.root,
        ...(opts.approver ? [] : [span('Approvers approve: approving meets the criterion', 'pp-muted')]),
        errors.root,
      ], 'pp-save-run'), 'link-document-form'),
      buttons: [{label: 'LINK', name: 'link-document-button', primary: true, enabled: valid,
        run: async () => {
          let target = choiceValue(folderFile) || path.value.value.trim();
          let folder: string | null = choiceValue(folderFile) ? program.space_id ?? null : null;
          const file = upload.files?.[0];
          if (file) {
            // uploaded once: a retry after a failed link links the file already in the folder
            if (uploaded?.file !== file)
              uploaded = {file, ...await uploadDocument(program, file.name, new Uint8Array(await file.arrayBuffer()))};
            target = uploaded.path;
            folder = uploaded.folder;
          }
          linked = await linkDocumentTo(program, {path: target, folder, criterionId: criterion?.id ?? null,
            kind: choiceValue(kind), title: title.value.value.trim() || title.placeholder,
            version: version.value.value.trim() || undefined, status: choiceValue(status) as DocumentStatus});
        },
        then: async () => {
          const r = linked;
          const doc = describeDocument({title: r.document.title, version: r.document.doc_version,
            status: r.document.status as DocumentStatus});
          const state = r.criterion ? `; ${r.criterion.code} ${r.criterion.state.replace('_', ' ')}` : '';
          grok.shell.info([`${program.code}: ${doc} linked${state}`, ...r.notes].join('\n'));
          await opts.onLinked?.(r);
        }}],
    };
  });
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
    const pre = named(document.createElement('pre'), 'document-text');
    pre.className = 'pp-card-text';
    pre.textContent = body;
    openDialog({title: describeDocument({title: d.title, version: d.doc_version, status: d.status as DocumentStatus}),
      name: 'document-viewer', width: 640}, () => ({content: divV([span(d.path, 'pp-muted'), pre])}));
    return;
  }
  DG.Utils.download(name, new Blob([(await grok.dapi.files.readAsBytes(d.path)).slice()]));
}
