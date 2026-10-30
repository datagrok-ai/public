# @datagrok-libraries/ems

UI and data layer for entity-mapped domain tables (EMS), built on `@datagrok-libraries/u2`.

```ts
import {domains} from '@datagrok-libraries/ems/src/dg';
import '@datagrok-libraries/ems/src/dg/styles.js';

const issues = await domains.table('grit.issue');
export function gritApp() { return issues.app(); }
```

## Layout

```
src/core/   platform-free: DomainSource, the backend seam (emsBackends), EditState,
            SharedSession, MemoryDomainBackend, the u2-domain-source tag
src/dg/     platform layer: the domains.* facade, the controls, DomainApp, authoring
css/        domain.css, binding-dialog.css, manifest-editor.css
docs/       dg.md, domain-backend-contract.md, recipes/
```

Dependencies run `ems → u2 → datagrok-api`; u2 never imports this library.

## Registration

`installEms()` adds the source and every `u2-domain-*` control to a registry (the global one by
default) — once, as a `registry.use()` plugin that also brings the ems unit of work to every spec
rendered over it — and answers the uninstall. The tags keep their `u2-` names so stored specs
still resolve.

## Development

```bash
pnpm -C libraries/ems run build
pnpm -C libraries/ems test
```

Browser checks reuse u2's harness (`../u2/e2e`) against a running local client:

```bash
pnpm -C libraries/ems run e2e:local          # the `domain` and `session` lanes over gallery/
pnpm -C libraries/ems run e2e:stand-domain   # the domain apps on a stand (U2_STAND_URL)
```

`gallery/index.html` hosts the two pages over u2's built `src/` through an import map; run
`pnpm test` in both libraries first, which compiles them into `.test-build/`, where `e2e/gallery.mjs`
serves the `.js` from.
