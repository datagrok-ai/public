/* What a domain control says when the source reports a failure. Over the platform backend the
   js-api editor never throws from `save()`: it balloons, maps the server's column errors onto
   the cells (which `errorOf` surfaces per field) and resolves a 409 through the platform's own
   reload/overwrite dialog — so nothing here duplicates that. What reaches `DomainSource.error`
   is the memory backend's refusal (the gallery, the tests) or a failed load, and both are one
   message; a 403 additionally drops the access caches so the affordances re-gate. */
import * as grok from 'datagrok-api/grok';
import {backends} from '../../sources/backends.js';
import {notify} from '../../components/display/notify.js';
import {DgDomainBackend} from './backend.js';

export class DomainErrors {
  static codeOf(e: unknown): string {
    const code = (e as {code?: unknown} | null)?.code;
    return typeof code === 'string' ? code : '';
  }

  static message(e: unknown): string {
    return e instanceof Error ? e.message : String((e as {message?: unknown} | null)?.message ?? e);
  }

  /** One balloon; after a 403 the affordances were built from an access snapshot the server just
   * contradicted, so every cache of it goes. */
  static report(e: unknown): void {
    if (DomainErrors.codeOf(e) === 'forbidden') {
      grok.dapi.domains.invalidateUiCaches();
      if (backends.domain instanceof DgDomainBackend)
        backends.domain.invalidate();
    }
    notify.error(DomainErrors.message(e));
  }
}
