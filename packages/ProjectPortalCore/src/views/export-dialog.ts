/** A loaded portfolio as its portfolio file: shown in a dialog, with a download. Opened from the
 * board's Export button; `exportPortfolio(code)` returns the same JSON. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {exportPortfolioConfig} from '../profile';

export async function openExportDialog(code: string): Promise<void> {
  const json = JSON.stringify(await exportPortfolioConfig(code), null, 2);
  const text = document.createElement('textarea');
  text.textContent = json;
  text.readOnly = true;
  text.style.cssText = 'width: 100%; height: 420px; font-family: monospace; font-size: 11px;';
  text.dataset.u2Name = 'export-json';
  const dlg = ui.dialog({title: `${code}: portfolio file`}).add(text);
  dlg.root.dataset.u2Name = 'export-dialog';
  dlg.addButton('Download', () => DG.Utils.download(`${code.toLowerCase()}.json`, json, 'application/json'));
  dlg.show({width: 640});
}
