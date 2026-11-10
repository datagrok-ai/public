/**
 * Regression guard for `mapping-pane-unit-dual-mode`.
 *
 * Each source unit (Time / Concentration / Dose / Body weight) carries TWO synchronised
 * pickers: a HORIZONTAL one inline on the value's row (`unitRow`), shown in the wide
 * (label-left) form layout, and a VERTICAL one as its own captioned row under the value,
 * shown in the narrow (label-on-top / `ui-form-condensed`) layout. `applyUnitPickerMode`
 * toggles which is visible off the form's `ui-form-condensed` class; `linkUnitPair` keeps
 * the pair equal so every consumer keeps reading the original (horizontal) input.
 *
 * This exists because a single inline picker orphans in the narrow layout (the reported
 * bug) — the two-picker crutch is the shipped fix. These tests pin (a) the visibility
 * toggle on the REAL pane in both layouts, (b) the pair starting equal, and (c) the sync
 * wiring both directions through the DG API.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {mountMappingPane} from '../views/nca-mapping-pane';
import {loadFixture} from './test-utils';

const TIMEOUT = 90000;
// Vertical pickers carry these captions; horizontal pickers carry an empty caption and
// sit inside a `.nca-unit-row` divH. Order matches the horizontal pickers' DOM order.
const V_CAPTIONS = ['Time unit', 'Concentration unit', 'Dose unit', 'Body weight unit'];

category('NCA: dual-mode unit pickers (mapping-pane-unit-dual-mode)', () => {
  test('real pane: horizontal shown wide, vertical shown narrow, pair equal', async () => {
    const df = await loadFixture('theoph.csv');
    const view = grok.shell.addTableView(df);
    try {
      mountMappingPane(view, df);
      await DG.delay(600);
      const paneRoot = document.querySelector('.nca-mapping-pane') as HTMLElement | null;
      expect(paneRoot !== null, true);

      const dockNode = paneRoot!.closest('.d4-dock-node') as HTMLElement | null;
      const setW = async (w: number): Promise<void> => {
        if (dockNode) dockNode.style.width = `${w}px`;
        if (paneRoot!.parentElement) paneRoot!.parentElement.style.width = `${w}px`;
        paneRoot!.style.width = `${w}px`;
        window.dispatchEvent(new Event('resize'));
        await DG.delay(500);
      };
      const vis = (el: HTMLElement): boolean => el.offsetParent !== null;
      const hPickers = (): HTMLElement[] => Array.from(
        paneRoot!.querySelectorAll('.nca-unit-row .nca-unit-picker')) as HTMLElement[];
      const vPickerFor = (caption: string): HTMLElement | undefined =>
        (Array.from(paneRoot!.querySelectorAll('.ui-input-root')) as HTMLElement[]).find((r) =>
          (r.querySelector('.ui-input-label')?.textContent ?? '').trim() === caption);
      const vPickers = (): HTMLElement[] =>
        V_CAPTIONS.map(vPickerFor).filter((r): r is HTMLElement => r != null);

      await setW(420); // WIDE / label-left
      expect(hPickers().filter(vis).length, 4); // all inline pickers shown
      expect(vPickers().filter(vis).length, 0); // vertical rows hidden

      // Pair starts equal: each horizontal picker's value equals its vertical twin's.
      const hs = hPickers();
      for (let i = 0; i < V_CAPTIONS.length; i++) {
        const hSel = hs[i]?.querySelector('select') as HTMLSelectElement | null;
        const vSel = vPickerFor(V_CAPTIONS[i])?.querySelector('select') as HTMLSelectElement | null;
        expect(hSel !== null && vSel !== null, true);
        expect(hSel!.value, vSel!.value);
      }

      await setW(150); // NARROW / label-on-top
      expect(hPickers().filter(vis).length, 0); // inline pickers hidden
      expect(vPickers().filter(vis).length, 4); // vertical rows shown
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});

  // The `linkUnitPair` wiring the pane uses, verified through the DG API (`.value =`
  // fires onChanged by contract). Raw DOM `<select>` edits do NOT drive DG's Dart-backed
  // value, so the pane's sync is verified here via the API, not by faking DOM events.
  test('linkUnitPair keeps the pair equal both directions (DG API)', async () => {
    const h = ui.input.choice<string>('', {value: 'h', items: ['h', 'min', 'day']});
    const v = ui.input.choice<string>('Time unit', {value: 'h', items: ['h', 'min', 'day']});
    h.onChanged.subscribe(() => {if (v.value !== h.value) v.value = h.value;});
    v.onChanged.subscribe(() => {if (h.value !== v.value) h.value = v.value;});
    h.value = 'min';
    await DG.delay(50);
    expect(v.value, 'min'); // H -> V
    v.value = 'day';
    await DG.delay(50);
    expect(h.value, 'day'); // V -> H
  }, {timeout: TIMEOUT});

  // Collapse-safety (decision-audit regression): a mode flip (resize-driven) must NOT
  // pop a vertical unit-picker row back into view under a COLLAPSED section header.
  // `applyUnitPickerMode` composes the mode intent with each section's collapse state
  // via the section `reapply` handle; without that fix, narrowing while "Profile data"
  // is collapsed un-hid the Time/Concentration vertical rows individually.
  test('collapse-safe: narrowing does not un-hide a collapsed section row', async () => {
    const df = await loadFixture('theoph.csv');
    const view = grok.shell.addTableView(df);
    try {
      mountMappingPane(view, df);
      await DG.delay(600);
      const paneRoot = document.querySelector('.nca-mapping-pane') as HTMLElement | null;
      expect(paneRoot !== null, true);
      const dockNode = paneRoot!.closest('.d4-dock-node') as HTMLElement | null;
      const setW = async (w: number): Promise<void> => {
        if (dockNode) dockNode.style.width = `${w}px`;
        if (paneRoot!.parentElement) paneRoot!.parentElement.style.width = `${w}px`;
        paneRoot!.style.width = `${w}px`;
        window.dispatchEvent(new Event('resize'));
        await DG.delay(400);
      };
      const rootByCaption = (caption: string): HTMLElement | undefined =>
        (Array.from(paneRoot!.querySelectorAll('.ui-input-root')) as HTMLElement[]).find((r) =>
          (r.querySelector('.ui-input-label')?.textContent ?? '').trim() === caption);

      await setW(420); // wide: vertical rows are the hidden ones
      const profileHeader = (Array.from(paneRoot!.querySelectorAll('h2')) as HTMLElement[])
        .find((h) => (h.textContent ?? '').includes('Profile data'));
      expect(profileHeader != null, true);
      profileHeader!.click(); // collapse "Profile data" (holds Time + Concentration pickers)
      await DG.delay(100);
      await setW(150); // narrow flip — the resize path that used to pop the rows back in

      // Time / Concentration vertical pickers are in the COLLAPSED section → stay hidden.
      expect(rootByCaption('Time unit')!.offsetParent === null, true);
      expect(rootByCaption('Concentration unit')!.offsetParent === null, true);
      // Control: Dose unit is in Context (still expanded) → visible in narrow mode.
      expect(rootByCaption('Dose unit')!.offsetParent !== null, true);
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});
});
