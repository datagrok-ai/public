/**
 * Moved into the platform. {@link DomainFrameEditor} and the editing-state
 * contract now ship inside `datagrok-api` (`src/ui/domains/domains-editor.ts`),
 * so every session carries them with no package installed — this module keeps
 * the import path, and the names below are the SAME classes `DG.*` carries,
 * never a second copy.
 *
 * The generic ones (the service-column names, the cell validator, the editor
 * rollup) are statics of `DG.DomainFrameEditor` on the platform surface; here
 * they keep the flat names every existing import already spells.
 *
 * @module frame-editor
 */

import {DomainFrameEditor} from 'datagrok-api/dg';

export {DomainFrameEditor, acquireDomainContext} from 'datagrok-api/dg';

export type {AnyDomainTableClient, DomainCellError, DomainFrameEditorOptions,
  DomainPendingOp, DomainRowState, DomainSaveResult, IDomainTableContext, IEditorHost} from 'datagrok-api/dg';

export const STATE_COLUMN = DomainFrameEditor.STATE_COLUMN;
export const CHANGES_COLUMN = DomainFrameEditor.CHANGES_COLUMN;
export const ERRORS_COLUMN = DomainFrameEditor.ERRORS_COLUMN;
export const SERVICE_COLUMNS = DomainFrameEditor.SERVICE_COLUMNS;
export const validateCellValue = DomainFrameEditor.validateCellValue;
export const isReferenceProperty = DomainFrameEditor.isReferenceProperty;
export const editorsOf = DomainFrameEditor.editorsOf;
