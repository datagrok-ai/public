/* eslint-disable max-len */
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {
  applyHistoryEntry,
  buildChemEnumPanel,
  buildHistoryEntry,
  ChemEnumDialogState,
  ChemEnumHistoryEntry,
} from './pt-chem-enum-dialog';

export const markushDefaultsPropName = 'MarkushDefaults';
export const monomersPathPropName = 'MonomersPath';

/** Serializes a panel state to the JSON string persisted in package settings — same shape as a history entry. */
export function serializeChemEnumState(state: ChemEnumDialogState): string {
  return JSON.stringify(buildHistoryEntry(state));
}

/** Parses a settings/defaults JSON blob into an entry; returns null on any parse or shape failure. */
export function parseChemEnumDefaults(json: string | null | undefined): ChemEnumHistoryEntry | null {
  if (!json) return null;
  try {
    const parsed = JSON.parse(json);
    if (!Array.isArray(parsed?.cores) || !parsed.cores.every((s: unknown) => typeof s === 'string') ||
      typeof parsed?.rGroups !== 'object' || parsed.rGroups === null)
      return null;
    return parsed as ChemEnumHistoryEntry;
  } catch {
    return null;
  }
}

export async function getMarkushDefaults(): Promise<ChemEnumHistoryEntry | null> {
  try {
    const ps = (_package.settings ?? await _package.getSettings()) as unknown as Map<string, any> | Record<string, any>;
    if (!ps)
      return null;
    // this hack is needed because api says it will return Map, but it actually returns object...
    const value = ps instanceof Map ? ps.get(markushDefaultsPropName) : ps[markushDefaultsPropName];
    return parseChemEnumDefaults(value ?? null);
  } catch (e) {
    console.error(e);
    return null;
  }
}

export async function markushSettingsEditorWidget(propList: DG.Property[]): Promise<DG.Widget> {
  const w = DG.Widget.fromRoot(ui.divV([], {style: {overflow: 'visible'}}));
  w.root.appendChild(ui.h1('Markush Enumerator Defaults', {style: {marginBottom: '6px', marginTop: '6px'}}));

  const defaultsProp = propList.find((p) => p.name === markushDefaultsPropName);
  if (defaultsProp) {
    const summary = ui.divText('');
    const updateSummary = () => {
      const curValue: string | null = defaultsProp.get(null); // null is passed here because there is an getter override in the core
      const entry = parseChemEnumDefaults(curValue);
      summary.textContent = entry ? `Configured: ${entry.summary || '(custom set)'}` : 'No defaults configured — the shipped CSV files are used.';
    };
    updateSummary();

    const editButton = ui.button('Edit defaults', async () => {
      let rdkit: RDModule;
      try {
        rdkit = await getRdKitModule();
      } catch (e) {
        grok.shell.error('Cannot open the editor — RDKit is unavailable. Check console for more details');
        console.error(e);
        return;
      }
      // The cores / R-groups editor is too rich for the narrow settings pane, so it is opened in a
      // dialog; OK stages the value via prop.set, the same way HitTriage's inputs set on change.
      const panel = buildChemEnumPanel(rdkit, null, 'dialog', {hideAppendToTable: true});
      const curValue: string | null = defaultsProp.get(null); // null is passed here because there is an getter override in the core
      const entry = parseChemEnumDefaults(curValue);
      if (entry) {
        applyHistoryEntry(entry, panel.state, rdkit);
        panel.refresh();
      }
      const dialog = ui.dialog({title: 'Edit Markush defaults'})
        .add(panel.root)
        .onOK(() => {
          try {
            defaultsProp.set(null, serializeChemEnumState(panel.state)); // FYI for future devs: null is passed here because there is an setter override in the core, the temp prop is set in backend and it is saved when a save button is clicked
            updateSummary();
          } catch (e) {
            grok.shell.error('Error saving Markush defaults. Check console for more details');
            console.error(e);
          }
        });
      panel.bindActionButton(dialog.getButton('OK') as HTMLButtonElement);
      dialog.show({resizable: true, width: 960});
    });

    const clearButton = ui.button('Clear defaults', () => {
      try {
        defaultsProp.set(null, ''); // FYI for future devs: null is passed here because there is an setter override in the core, the temp prop is set in backend and it is saved when a save button is clicked
        updateSummary();
      } catch (e) {
        grok.shell.error('Error clearing Markush defaults. Check console for more details');
        console.error(e);
      }
    });
    ui.tooltip.bind(clearButton, 'Remove the configured defaults; the app falls back to the shipped CSV files');

    w.root.appendChild(summary);
    w.root.appendChild(ui.divH([editButton, clearButton], {style: {gap: '8px', marginTop: '6px'}}));
  }

  const monomersProp = propList.find((p) => p.name === monomersPathPropName);
  if (monomersProp) {
    const defaultTooltip = monomersProp.description ?? 'Path to additional monomer libraries';
    let curTooltip = defaultTooltip;
    const curValue: string = monomersProp.get(null) ?? monomersProp.defaultValue; // null is passed here because there is an getter override in the core
    const monomersHeader = ui.h1('Monomers', {style: {marginBottom: '6px', marginTop: '6px'}});
    const monomersInput = ui.input.string('Monomers path', {value: curValue});
    ui.tooltip.bind(monomersInput.input, () => curTooltip);
    DG.debounce(monomersInput.onChanged, 200).subscribe((_) => {
      const value = monomersInput.value;
      try {
        if (!value) {
          monomersProp.set(null, monomersProp.defaultValue);
          curTooltip = defaultTooltip + '<br><br> <b style="color: red">Cannot be empty.</b>';
          monomersInput.input.classList.add('d4-invalid');
          return;
        }
        monomersProp.set(null, value);
        curTooltip = defaultTooltip;
        monomersInput.input.classList.remove('d4-invalid');
      } catch (e) {
        grok.shell.error('Error Setting monomers path. Check console for more details');
        console.error(e);
      }
    });
    w.root.appendChild(monomersHeader);
    w.root.appendChild(monomersInput.root);
  }

  setTimeout(() => w.root.parentElement?.style && (w.root.parentElement.style.overflow = 'visible'), 200);
  return w;
}
