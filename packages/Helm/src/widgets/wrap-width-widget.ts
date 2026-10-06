import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';

// The wrapping width is a process-wide setting inside hwe, read at layout time
// by every HelmService / editor / adapter (see @datagrok-libraries/hwe API.md
// §4.3). That is what lets one input here change both the grid cells and the
// editor without hunting down the several service instances the package owns.
import {getWrapSettings, setMonomersPerRow, onWrapSettingsChanged} from '@datagrok-libraries/hwe';

import {_package} from '../package';

/** Widest row the input allows. Well past any readable structure. */
const MAX_MONOMERS_PER_ROW = 100;

/** Repaint every open grid so HELM cells re-render at the new width. */
function invalidateHelmGrids(): void {
  // hwe's own layout / render caches key on a setting version and have already
  // missed by the time we get here, so a plain repaint is enough.
  for (const view of grok.shell.tableViews) {
    try {
      view.grid?.invalidate();
    } catch (_err) {
      // A view can be closing while we iterate; a failed repaint is not worth
      // failing the whole widget over.
    }
  }
}

/**
 * Session-scoped control over where HELM structures break into a new line.
 *
 * Deliberately does NOT write the package property: the `MonomersPerRow`
 * package setting stays the default that every session starts from, and this
 * widget adjusts the current session only. "Reset" puts the session back onto
 * that default rather than clearing it.
 *
 * @return {DG.Widget} The row-width editor.
 */
export function getWrapWidthWidget(): DG.Widget {
  const input = ui.input.int('Monomers per row', {
    value: getWrapSettings().monomersPerRow,
    min: 1,
    max: MAX_MONOMERS_PER_ROW,
    tooltipText: 'Number of monomers drawn per row before the structure wraps. ' +
      'Counts nucleotides for RNA/DNA. Applies to this session only. The default comes from the Helm package property',
    onValueChanged: (value: number | null) => {
      if (value === null || !Number.isFinite(value)) return;
      setMonomersPerRow(value);
      invalidateHelmGrids();
    },
  });

  const packageDefault = _package.defaultMonomersPerRow;
  const resetLink = ui.link('Reset to package default', () => {
    setMonomersPerRow(packageDefault);
    invalidateHelmGrids();
  }, `Back to ${packageDefault} — the MonomersPerRow package property`, {style: {marginTop: '4px'}});

  const widget = new DG.Widget(ui.divV([input.root, resetLink]));

  // Keep the input honest if the setting is changed from anywhere else (the
  // editor, another widget, a script). `setMonomersPerRow` clamps and rounds,
  // so this is also how the input reflects a value it did not produce itself.
  const off = onWrapSettingsChanged((settings) => {
    if (input.value !== settings.monomersPerRow) input.value = settings.monomersPerRow;
  });
  widget.sub(new Subscription(off));

  return widget;
}
