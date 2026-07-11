import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parseProjectVocabulary} from './publish-settings';

/**
 * Custom editor for the Proteomics package settings (Plugins → Proteomics →
 * Settings). Replaces the auto-generated form so `projectVocabulary` gets a
 * proper add / edit / remove list instead of a single comma-separated text box.
 *
 * Follows the platform contract used by HitTriage / SequenceTranslator: read the
 * live value with `prop.get(null)` and stage edits with `prop.set(null, value)`;
 * the Settings pane's SAVE button commits the staged values. (Direct
 * `Package.setSettings()` is NOT the write path here — the property staging is.)
 *
 * Because a custom editor replaces the ENTIRE settings form, every Proteomics
 * setting is rendered here, not just the vocabulary — otherwise the others would
 * vanish from the pane.
 */
export function proteomicsSettingsEditorWidget(propList: DG.Property[]): DG.Widget {
  const findProp = (name: string): DG.Property | undefined => propList.find((p) => p.name === name);

  const inputs: DG.InputBase[] = [];

  const addStringInput = (name: string, caption: string): void => {
    const prop = findProp(name);
    if (!prop) return;
    const cur = (prop.get(null) as string | null) ?? '';
    const input = ui.input.string(caption, {value: cur});
    input.onChanged.subscribe(() => {
      try { prop.set(null, input.value ?? ''); } catch (e) { console.error(e); }
    });
    inputs.push(input);
  };

  const addBoolInput = (name: string, caption: string): void => {
    const prop = findProp(name);
    if (!prop) return;
    const cur = prop.get(null);
    const input = ui.input.bool(caption, {value: cur === true || cur === 'true'});
    input.onChanged.subscribe(() => {
      try { prop.set(null, !!input.value); } catch (e) { console.error(e); }
    });
    inputs.push(input);
  };

  addStringInput('reviewSpaceName', 'Review Space Name');
  addStringInput('reviewNamePrefix', 'Review Name Prefix');
  addBoolInput('verifyPublishedDashboard', 'Verify Published Dashboard');

  const root = ui.divV([ui.form(inputs)]);

  // ── Project vocabulary — add / edit / remove list ───────────────────────────
  const vocabProp = findProp('projectVocabulary');
  if (vocabProp) {
    // Serialize on every edit; parseProjectVocabulary trims + dedupes on read.
    let items = parseProjectVocabulary(vocabProp.get(null));
    const commit = (): void => {
      try { vocabProp.set(null, items.join(', ')); } catch (e) { console.error(e); }
    };

    const rowsHost = ui.divV([]);
    const renderRows = (): void => {
      ui.empty(rowsHost);
      if (items.length === 0) {
        rowsHost.appendChild(ui.divText('No projects yet — add one below.',
          {style: {color: 'var(--grey-4)', fontStyle: 'italic', margin: '4px 0'}}));
      }
      items.forEach((val, i) => {
        const input = ui.input.string('', {value: val});
        input.input.style.width = '260px';
        input.onChanged.subscribe(() => { items[i] = (input.value ?? '').trim(); commit(); });
        const del = ui.icons.delete(() => { items.splice(i, 1); renderRows(); commit(); }, 'Remove project');
        rowsHost.appendChild(ui.divH([input.input, del], {style: {alignItems: 'center', gap: '8px', margin: '2px 0'}}));
      });
    };

    const addBtn = ui.button('+ Add project', () => { items.push(''); renderRows(); });
    renderRows();

    root.appendChild(ui.h1('Projects', {style: {marginTop: '12px', marginBottom: '4px'}}));
    root.appendChild(ui.divText('Analysts share each analysis under one of these. Add, rename, or remove names here.',
      {style: {color: 'var(--grey-5)', marginBottom: '6px'}}));
    root.appendChild(rowsHost);
    root.appendChild(addBtn);
  }

  return DG.Widget.fromRoot(root);
}
