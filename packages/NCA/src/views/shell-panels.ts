/**
 * The platform side panel around an analysis view (user-directed 2026-10-12): the Browse tree in
 * front, the toolbox closed — and the custom event a docked compute fires so a Browse-tree preview
 * can take keyboard focus back. A LEAF module (grok only): the open paths (nca-launch,
 * nca-analysis-open), the mapping pane and the Browse tree all import it without an import cycle
 * and without the tree having to load the mapping pane for a string constant.
 */

import * as grok from 'datagrok-api/grok';

/** Custom-event id fired (with the analysis view) once a compute's results are docked. The
 *  Browse tree listens so a tree preview can take keyboard focus back from the docked grid. */
export const NCA_RESULTS_PRESENTED_EVENT = 'nca-results-presented';

/** The DiffStudio pair — toolbox off, the Browse tree in front — plus the platform's auto-show off,
 *  so opening a view does not bring the toolbox back (user-directed 2026-10-12; a GLOBAL setting,
 *  deliberately not restored). Browse and the toolbox share the side panel, and a table view docks
 *  its toolbox as a tab OVER Browse, so hiding the toolbox alone left that tab in front.
 *  Idempotent. Called as the FIRST step of opening an analysis (before the study and the table
 *  load), right before the view is added, and again once a loaded analysis has docked its results
 *  (the docking brought the toolbox back on a large table — QA 2026-10-12). */
export function keepBrowseInFront(): void {
  grok.shell.windows.autoShowToolbox = false;
  grok.shell.windows.showToolbox = false;
  grok.shell.windows.showBrowse = true;
}
