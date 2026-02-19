/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';

const MAX_SUGGESTIONS = 20;

/** Shows a dialog for selecting monomers with autocomplete and tag-based display.
 * @returns comma-separated monomer symbols, or null if cancelled */
export async function showMonomerSelectionDialog(
  monomerLib: IMonomerLib, polymerType: PolymerType, presetMonomers?: string[],
): Promise<string[] | null> {
  return new Promise<string[] | null>((resolve) => {
    const helmType: HelmType = polymerTypeToHelmType(polymerType);
    const allSymbols = monomerLib.getMonomerSymbolsByType(polymerType);

    const selectedMonomers: string[] = presetMonomers ? [...presetMonomers] : [];

    const tagsHost = ui.div([], {style: {display: 'flex', flexWrap: 'wrap', gap: '4px', marginTop: '8px', maxWidth: '400px'}});
    const input = ui.input.string('Monomers', {value: ''});
    const inputEl = input.input as HTMLInputElement;
    inputEl.setAttribute('autocomplete', 'off');
    inputEl.placeholder = 'Type to search...';

    let currentMenu: DG.Menu | null = null;
    let menuItems: HTMLElement[] = [];
    let highlightedIndex = -1;

    function renderTags(): void {
      tagsHost.innerHTML = '';
      for (const symbol of selectedMonomers) {
        const removeBtn = ui.iconFA('times', () => {
          const idx = selectedMonomers.indexOf(symbol);
          if (idx >= 0) {
            selectedMonomers.splice(idx, 1);
            renderTags();
          }
        });
        removeBtn.style.marginLeft = '4px';
        removeBtn.style.cursor = 'pointer';

        const tag = ui.div([ui.span([symbol]), removeBtn], {
          style: {
            display: 'inline-flex', alignItems: 'center',
            padding: '2px 6px', borderRadius: '4px',
            backgroundColor: 'var(--grey-2)', border: '1px solid var(--grey-3)',
            fontSize: '12px', cursor: 'default',
          },
        });
        // Tooltip on hover
        tag.addEventListener('mouseenter', (e) => {
          const tooltip = monomerLib.getTooltip(helmType, symbol);
          ui.tooltip.show(tooltip, tag.getBoundingClientRect().left, tag.getBoundingClientRect().bottom + 16);
        });
        tag.addEventListener('mouseleave', () => { ui.tooltip.hide(); });

        tagsHost.appendChild(tag);
      }
    }

    function addMonomer(symbol: string): void {
      if (!selectedMonomers.includes(symbol)) {
        selectedMonomers.push(symbol);
        renderTags();
      }
      inputEl.value = '';
      hideMenu();
      inputEl.focus();
    }

    function hideMenu(): void {
      if (currentMenu) {
        currentMenu.hide();
        currentMenu = null;
      }
      menuItems = [];
      highlightedIndex = -1;
    }

    function getSuggestions(query: string): {symbol: string, monomer: Monomer | null}[] {
      const q = query.toLowerCase();
      const results: {symbol: string, monomer: Monomer | null, rank: number}[] = [];

      for (const symbol of allSymbols) {
        if (selectedMonomers.includes(symbol))
          continue;
        const monomer = monomerLib.getMonomer(polymerType, symbol);
        const symLower = symbol.toLowerCase();
        const nameLower = monomer?.name?.toLowerCase() ?? '';

        if (symLower.startsWith(q))
          results.push({symbol, monomer, rank: 0});
        else if (symLower.includes(q))
          results.push({symbol, monomer, rank: 1});
        else if (nameLower.includes(q))
          results.push({symbol, monomer, rank: 2});
      }

      results.sort((a, b) => a.rank - b.rank || a.symbol.localeCompare(b.symbol));
      return results.slice(0, MAX_SUGGESTIONS);
    }

    function showSuggestions(): void {
      hideMenu();
      const query = inputEl.value.trim();
      if (!query)
        return;

      const suggestions = getSuggestions(query);
      if (suggestions.length === 0)
        return;

      currentMenu = DG.Menu.popup();
      const maxElement = suggestions.reduce((max, s) => {
        const label = s.monomer?.name ? `${s.symbol} - ${s.monomer.name}` : s.symbol;
        if (max.length < label.length)
          return label;
        return max;
      }, '');
      currentMenu.item(maxElement, () => {}); // Dummy item to set menu width

      const causedBy = new MouseEvent('mousemove', {clientX: inputEl.getBoundingClientRect().left, clientY: inputEl.getBoundingClientRect().bottom});
      currentMenu.show({causedBy: causedBy,
        y: inputEl.offsetHeight + inputEl.offsetTop, x: inputEl.offsetLeft});

      // collect menu items for keyboard navigation
      setTimeout(() => {
        currentMenu?.clear();
        for (const s of suggestions) {
          const label = s.monomer?.name ? `${s.symbol} - ${s.monomer.name}` : s.symbol;
          currentMenu?.item(label, () => { addMonomer(s.symbol); });
        }

        const menuRoot = document.querySelector('.d4-menu-popup:last-of-type');
        if (menuRoot) {
          menuItems = Array.from(menuRoot.querySelectorAll('.d4-menu-item')) as HTMLElement[];
          menuItems.forEach((item) => {
            item.addEventListener('mouseenter', () => {
              const symbol = item.textContent?.split(' - ')[0] ?? '';
              if (!symbol)
                return;
              const tooltip = monomerLib.getTooltip(helmType, symbol);
              ui.tooltip.show(tooltip, item.getBoundingClientRect().right + 10, item.getBoundingClientRect().bottom + 16);
            });
          });
        }
        highlightedIndex = -1;
      }, 0);
    }

    function updateHighlight(newIndex: number): void {
      if (menuItems.length === 0)
        return;
      if (highlightedIndex >= 0 && highlightedIndex < menuItems.length)
        menuItems[highlightedIndex].classList.remove('d4-menu-item-hover');
      highlightedIndex = newIndex;
      if (highlightedIndex >= 0 && highlightedIndex < menuItems.length) {
        menuItems[highlightedIndex].classList.add('d4-menu-item-hover');
        menuItems[highlightedIndex].scrollIntoView({block: 'nearest'});
      }
    }

    inputEl.addEventListener('input', () => { showSuggestions(); });

    inputEl.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        if (menuItems.length > 0) {
          const next = (highlightedIndex + 1) % menuItems.length;
          updateHighlight(next);
        }
      } else if (e.key === 'ArrowUp') {
        e.preventDefault();
        if (menuItems.length > 0) {
          const prev = (highlightedIndex - 1 + menuItems.length) % menuItems.length;
          updateHighlight(prev);
        }
      } else if (e.key === 'Enter') {
        e.preventDefault();
        e.stopPropagation();
        if (highlightedIndex >= 0 && highlightedIndex < menuItems.length) {
          menuItems[highlightedIndex].click();
        } else {
          // If input exactly matches a symbol, add it directly
          const val = inputEl.value.trim();
          if (val && allSymbols.includes(val))
            addMonomer(val);
        }
      } else if (e.key === 'Escape') {
        hideMenu();
      }
    });

    renderTags();

    const dlg = ui.dialog({title: 'Select Monomers', showFooter: true})
      .add(ui.div([input.root, tagsHost], {style: {minWidth: '350px', minHeight: '200px'}}))
      .onOK(() => { resolve(selectedMonomers); })
      .onCancel(() => { resolve(null); })
      .show({resizable: true});
    // dlg.root.addEventListener('close', () => { hideMenu(); });

    inputEl.focus();
    setTimeout(() => { showSuggestions(); }, 0);
  });
}
