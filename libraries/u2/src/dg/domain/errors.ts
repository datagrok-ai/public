/* What a domain control says when the source reports a failure. Over the platform backend the
   js-api editor never throws from `save()`: it balloons, maps the server's column errors onto
   the cells (which `errorOf` surfaces per field) and resolves a 409 through the platform's own
   reload/overwrite dialog — so nothing here duplicates that. What reaches `DomainSource.error`
   is the memory backend's refusal (the gallery, the tests) or a failed load, and both are one
   message; a 403 additionally drops the access caches so the affordances re-gate. */
import * as grok from 'datagrok-api/grok';
import {backends} from '../../sources/backends.js';
import {Rows} from '../../sources/rows-like.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {notify} from '../../components/display/notify.js';
import {DgDomainBackend} from './backend.js';

/** How both backends word the veto on restoring a child whose parent is still deleted
 * (`repository.dart:2386`, `memory-domain.ts:297`): a column and the table it points at. Matched
 * by its shape wherever it stands, whatever the code says: over the platform backend the server's
 * `restrict` reaches this as the session's refusal text ("Cannot save: …"), whose code is
 * `'refused'` — the sentence is the only thing both paths share. */
const RESTRICT_PARENT = /Column "(.+?)" references a deleted row in "(.+?)"/;

export class DomainErrors {
  /** The seam's one code, from either backend: the memory one throws `DomainBackendError`, the
   * platform's `DomainError.code` is the server's `body.error` — so `'unsupported'` reads the
   * same whether the table refused the op here or the engine refused it there. */
  static codeOf(e: unknown): string {
    const code = (e as {code?: unknown} | null)?.code;
    return typeof code === 'string' ? code : '';
  }

  /** `source` is what the sentence is worded from where the refusal is about one of its rows. */
  static message(e: unknown, source?: DomainSource): string {
    const raw = e instanceof Error ? e.message : String((e as {message?: unknown} | null)?.message ?? e);
    const parts = RESTRICT_PARENT.exec(raw);
    if (parts === null)
      return raw;
    // the backend names a column and a table; the user asked to bring a ROW back, and the way
    // forward is the parent, not the column
    const prop = source?.schema.properties.find((p) => p.name === parts[1]);
    const target = (prop?.friendlyName ?? parts[2]).replace(/_/g, ' ').toLowerCase();
    const staged = (source?.pending() ?? []).filter((r) => r[Rows.STATE] === 'restored');
    const name = source?.schema.info.nameColumn;
    const row = staged.length !== 1 ? null : (name == null ? null : staged[0][name]) ?? staged[0].id;
    return `Cannot restore${row === null ? '' : ` "${String(row)}"`}: its ${target} was deleted too. ` +
      `Restore the ${target} first.`;
  }

  /** One balloon; after a 403 the affordances were built from an access snapshot the server just
   * contradicted, so every cache of it goes. */
  static report(e: unknown, source?: DomainSource): void {
    if (DomainErrors.codeOf(e) === 'forbidden') {
      grok.dapi.domains.invalidateUiCaches();
      if (backends.domain instanceof DgDomainBackend)
        backends.domain.invalidate();
    }
    notify.error(DomainErrors.message(e, source));
  }
}
