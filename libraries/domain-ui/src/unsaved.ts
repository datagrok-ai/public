/**
 * Moved into the platform. The save / discard / cancel prompt now ships inside
 * `datagrok-api` (`src/ui/domains/domains-editor.ts`, next to the editor whose
 * batch it stands in front of) — this module keeps the import path.
 *
 * @module unsaved
 */

export {confirmDiscardChanges, promptUnsavedChanges} from 'datagrok-api/dg';

export type {UnsavedOutcome, UnsavedPromptOptions} from 'datagrok-api/dg';
