/**
 * The Reactions help page (without it, the help pane said "Help unavailable" for the app's dialogs). Datagrok's own
 * help has no page on reactions, and its help tree is not in this repository, so Chem ships one, help/reactions.md, as
 * other packages ship theirs (HitTriage, ClinicalCase). The help pane reads it from the package. The app's view and its
 * Reset to example, the Reaction panel's full screen, the Export panel and the reaction cell editor each point at their
 * section.
 */
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';

/** The page, in the package's files. */
export const REACTIONS_HELP_FILE = 'help/reactions.md';

/** The page's sections, as its headings' anchors. */
export type ReactionsHelpSection = 'the-reaction-sketcher' | 'the-reaction-panel' | 'editing-a-reaction-cell' |
  'export';

/** Chem's web root, ending in a slash. In the package-test bundle (a second bundle) package.ts's `_package` was never
 * given to the platform, so its `webRoot` throws: the published Chem's is read from one of its functions there. */
export function chemWebRoot(): string {
  let root = '';
  try {
    root = _package.webRoot;
  } catch {
    root = DG.Func.find({package: 'Chem', name: 'getRdKitModule'})[0]?.package?.webRoot ?? '';
  }
  return root === '' || root.endsWith('/') ? root : `${root}/`;
}

/** The page's address, at `section`. */
export function reactionsHelpUrl(section?: ReactionsHelpSection): string {
  return `${chemWebRoot()}${REACTIONS_HELP_FILE}${section ? `#${section}` : ''}`;
}
