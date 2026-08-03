/**
 * `@datagrok-libraries/domain-ui` — composable UI over entity-mapped domain
 * tables (`grok.dapi.domains`).
 *
 * Everything here is reflective: components take columns, labels, choices,
 * validation and permissions from the runtime domain registry, so they work on
 * ANY registered table without codegen.
 *
 * ```ts
 * import {DomainGrid} from '@datagrok-libraries/domain-ui';
 * const grid = await DomainGrid.create(grok.dapi.domains.table('grit.issue'));
 * grok.shell.newView('Issues', [grid.root]);
 * ```
 *
 * A whole browse/CRUD app is the same length — {@link DomainAppView} is the list
 * page (search, cards/grid, New, deep-linkable query) and opens
 * {@link DomainEntityAppView} for a row (form, detail tabs, history):
 *
 * ```ts
 * grok.shell.addView(new DomainAppView(grok.dapi.domains.table('grit.issue')));
 * ```
 *
 * @module domain-ui
 */

export * from './frame-editor';
export * from './domain-grid';
export * from './entity-list';
export * from './app-view';
export * from './handler';
export * from './unsaved';
export * from './styles';

// The core domain surface, re-exported so apps import from ONE place. These are
// the SAME classes `DG.*` carries (they resolve to the `datagrok-api/dg`
// external at bundle time) — this is a convenience alias, never a second copy.
export {
  DomainObjectHandler,
  DomainView,
  DomainQuery,
  DomainError,
  DomainValidationError,
  DomainVersionConflictError,
  DomainRestrictError,
  DomainForbiddenError,
  DomainNotFoundError,
  DomainRow,
  retryOnVersionConflict,
  splitDomainTable,
  cond, and, or,
} from 'datagrok-api/dg';

export type {
  DomainAction,
  DomainDetailTab,
  DomainViewOptions,
  DomainTableCapabilities,
  DomainTableInfo,
  DomainChildTableRef,
  DomainQuerySpec,
  DomainQueryParams,
  DomainTransactionOp,
  DomainColumnError,
  DomainRowErrors,
} from 'datagrok-api/dg';
