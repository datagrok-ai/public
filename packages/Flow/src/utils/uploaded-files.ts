/** Local-file upload support for Flow.
 *
 * Bytes dropped onto the canvas are held in an in-memory registry under a
 * temporary `pending:` id until the flow is saved. Saving persists each file
 * to the server's GUID-addressed blob store (a connectionless `DG.FileInfo` —
 * bytes go to `files/data/{id}`, the same mechanism Compute uses for
 * replayable file inputs of historical runs) and rewrites the node's `fileId`
 * to the real entity id. `readUploadedFile` checks the registry first, so an
 * unsaved flow replays from memory with no server round-trip. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {parseFlowBody} from '../serialization/flow-script-format';

export const MAX_UPLOAD_SIZE = 100 * 1024 * 1024; // 100 MB

const PENDING_PREFIX = 'pending:';
const pending = new Map<string, {name: string; bytes: Uint8Array}>();

export function isPendingFileId(fileId: string): boolean {
  return fileId.startsWith(PENDING_PREFIX);
}

/** Registers dropped bytes; returns the temporary id to store in the node.
 *  Throws when the file exceeds {@link MAX_UPLOAD_SIZE}. */
export function addPendingFile(name: string, bytes: Uint8Array): string {
  if (bytes.length > MAX_UPLOAD_SIZE) {
    throw new Error(`"${name}" is ${formatSize(bytes.length)} — ` +
      `uploaded files are limited to ${formatSize(MAX_UPLOAD_SIZE)}`);
  }
  const id = PENDING_PREFIX + newId();
  pending.set(id, {name, bytes});
  return id;
}

export function getPendingFile(fileId: string): {name: string; bytes: Uint8Array} | null {
  return pending.get(fileId) ?? null;
}

export function removePendingFile(fileId: string): void {
  pending.delete(fileId);
}

/** Uploads a pending file to the server blob store and saves its FileInfo
 *  entity (which makes it addressable and shareable); returns the saved
 *  FileInfo whose `id` replaces the pending id in the node. */
export async function persistPendingFile(fileId: string): Promise<DG.FileInfo> {
  const p = pending.get(fileId);
  if (!p)
    throw new Error(`No pending upload for "${fileId}"`);
  const fi = DG.FileInfo.fromBytes(p.name, p.bytes);
  try {
    await grok.dapi.files.write(fi); // connectionless FileInfo → GUID-addressed blob
  } catch (e) {
    throw new Error(`Uploading "${p.name}" failed: ${errMsg(e)}`);
  }
  let saved: DG.FileInfo;
  try {
    saved = await fi.save();
  } catch (e) {
    // Servers without the connection-less FileInfo fix (files_service.dart
    // dereferenced f.connection unconditionally) reject the entity
    // registration. The blob itself is stored and readable by GUID, so
    // degrade: replay and sharing-by-id work, permission sync is skipped.
    console.warn(`FuncFlow: could not register "${p.name}" as a FileInfo entity — ` +
      `uploaded-file permissions cannot be synced on this server. ${errMsg(e)}`);
    pending.delete(fileId);
    return fi;
  }
  if (saved?.id == null)
    throw new Error(`Registering "${p.name}" on the server returned no id`);
  pending.delete(fileId);
  return saved;
}

/** The persisted (non-pending) uploaded-file ids referenced by a flow body. */
export function uploadedFileIdsFromFlowBody(body: string): string[] {
  try {
    const {doc} = parseFlowBody(body);
    return (doc.nodes ?? [])
      .filter((n) => String(n.typeName ?? '').toLowerCase().includes('readuploadedfile'))
      .map((n) => n.inputValues?.['fileId'])
      .filter((id): id is string => typeof id === 'string' && id !== '' && !isPendingFileId(id));
  } catch {
    return [];
  }
}

/** Grants read on every uploaded-file blob referenced by the flow script to
 *  the groups the script is visible to — a shared flow must be able to read
 *  its files. Driven by the package's autostart share listener (so it works
 *  with no Flow view open) and called after every save from the editor. */
export async function syncFlowFilePermissions(script: DG.Script): Promise<void> {
  const fileIds = uploadedFileIdsFromFlowBody(script.script ?? '');
  if (fileIds.length === 0) return;
  try {
    const perms = await grok.dapi.permissions.get(script) as
      unknown as {view?: DG.Group[]; edit?: DG.Group[]};
    const groups = [...(perms.view ?? []), ...(perms.edit ?? [])];
    if (groups.length === 0) return;
    for (const id of fileIds) {
      // The sharer may lack rights on a file uploaded by someone else —
      // skip that file, keep granting the rest.
      const fi = await grok.dapi.entities.find(id).catch(() => null);
      if (fi == null) continue;
      for (const g of groups) {
        try {
          await grok.dapi.permissions.grant(fi, g, false);
        } catch (e) {
          console.warn(`FuncFlow: could not grant file "${id}" to "${g.friendlyName}": ${errMsg(e)}`);
        }
      }
    }
  } catch (e) {
    console.warn('FuncFlow: failed to sync uploaded-file permissions', e);
  }
}

function errMsg(e: any): string {
  return String(e?.message ?? e?.innerMessage ?? e);
}

/** The bytes behind a node's fileId: pending registry first, then the server
 *  blob store. Errors are rephrased into something the flow user can act on. */
export async function readUploadedFileBytes(fileId: string, fileName: string): Promise<Uint8Array> {
  const p = pending.get(fileId);
  if (p)
    return p.bytes;
  if (isPendingFileId(fileId)) {
    throw new Error(`Uploaded file "${fileName}" is no longer available in this session — ` +
      'drop it onto the canvas again');
  }
  try {
    return await grok.dapi.files.readAsBytes(fileId);
  } catch (e) {
    throw new Error(`Cannot read uploaded file "${fileName}" — it may have been deleted, ` +
      `or you don't have access to it (ask the flow's owner to re-share). ${errMsg(e)}`);
  }
}

const TEXT_TABLE_EXTS = new Set(['csv', 'tsv', 'txt']);

/** Every file extension {@link parseFileToDataFrame} can turn into a table:
 *  the CSV family and `.d42` natively, plus whatever `file-handler` functions
 *  are registered (xlsx, sdf, …). Drives the upload picker's `accept` filter. */
export function supportedUploadExtensions(): string[] {
  const exts = new Set([...TEXT_TABLE_EXTS, 'd42']);
  for (const f of DG.Func.find({tags: ['file-handler']})) {
    for (const e of String(f.options['ext'] ?? '').split(','))
      if (e.trim()) exts.add(e.trim().toLowerCase());
  }
  return [...exts].sort();
}

/** Parses raw file bytes into a DataFrame: CSV-family natively, `.d42`
 *  directly, anything else through the platform's registered `file-handler`
 *  functions (xlsx, sdf, ...). Multi-table files yield their first table. */
export async function parseFileToDataFrame(fileName: string, bytes: Uint8Array): Promise<DG.DataFrame> {
  const ext = (fileName.split('.').pop() ?? '').toLowerCase();
  let df: DG.DataFrame;
  if (TEXT_TABLE_EXTS.has(ext))
    df = DG.DataFrame.fromCsv(new TextDecoder().decode(bytes));
  else if (ext === 'd42')
    df = DG.DataFrame.fromByteArray(bytes);
  else {
    const handler = findFileHandler(ext);
    if (!handler)
      throw new Error(`No importer is registered for ".${ext}" files`);
    const input = handler.inputs[0];
    const isText = input?.propertyType === DG.TYPE.STRING;
    const arg = isText ? new TextDecoder().decode(bytes) : bytes;
    const tables: DG.DataFrame[] = await handler.apply({[input.name]: arg});
    if (tables == null || tables.length === 0)
      throw new Error(`"${fileName}" contains no tables`);
    df = tables[0];
  }
  if (df.name == null || df.name === '')
    df.name = fileName.replace(/\.[^.]+$/, '');
  return df;
}

function findFileHandler(ext: string): DG.Func | null {
  for (const f of DG.Func.find({tags: ['file-handler']})) {
    const exts = String(f.options['ext'] ?? '').split(',').map((s) => s.trim().toLowerCase());
    if (exts.includes(ext))
      return f;
  }
  return null;
}

function newId(): string {
  const c = globalThis.crypto as Crypto | undefined;
  return c?.randomUUID != null ? c.randomUUID() :
    `${Date.now().toString(36)}-${Math.random().toString(36).slice(2)}`;
}

function formatSize(bytes: number): string {
  return bytes >= 1024 * 1024 ? `${Math.round(bytes / (1024 * 1024))} MB` : `${Math.round(bytes / 1024)} KB`;
}
