/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLibHelper, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';

import {parseMonomerSymbolList} from './pt-placeholders-input';

const MAX_SUGGESTIONS = 20;

/** Shows a dialog for selecting monomers with autocomplete, tag-based display,
 * and bulk-add from monomer libraries or collections.
 * @returns monomer symbols array, or null if cancelled */
export async function showMonomerSelectionDialog(
  monomerLib: IMonomerLib, polymerType: PolymerType, presetMonomers?: string[],
  libHelper?: IMonomerLibHelper,
): Promise<string[] | null> {
  return new Promise<string[] | null>((resolve) => {
    const helmType: HelmType = polymerTypeToHelmType(polymerType);
    const allSymbols = monomerLib.getMonomerSymbolsByType(polymerType);

    const selectedMonomers: string[] = presetMonomers ? [...presetMonomers] : [];

    // --- Tags display (selected monomers) ---
    const tagsHost = ui.div([], {style: {display: 'flex', flexWrap: 'wrap', gap: '4px', marginTop: '4px', maxHeight: '150px', overflowY: 'auto'}});
    const countLabel = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-4)', marginTop: '4px'}});

    function updateCountLabel(): void {
      countLabel.textContent = selectedMonomers.length > 0 ? `${selectedMonomers.length} monomer(s) selected` : '';
    }

    // --- Autocomplete input ---
    const input = ui.input.string('Search', {value: ''});
    const inputEl = input.input as HTMLInputElement;
    inputEl.setAttribute('autocomplete', 'off');
    inputEl.placeholder = 'Type to search monomers...';

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
            updateCountLabel();
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
        tag.addEventListener('mouseenter', () => {
          const tooltip = monomerLib.getTooltip(helmType, symbol);
          ui.tooltip.show(tooltip, tag.getBoundingClientRect().left, tag.getBoundingClientRect().bottom + 16);
        });
        tag.addEventListener('mouseleave', () => { ui.tooltip.hide(); });

        tagsHost.appendChild(tag);
      }
      updateCountLabel();
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

    function addMonomers(symbols: string[]): void {
      let added = false;
      for (const symbol of symbols) {
        if (!selectedMonomers.includes(symbol)) {
          selectedMonomers.push(symbol);
          added = true;
        }
      }
      if (added)
        renderTags();
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
        if (menuRoot)
          menuItems = Array.from(menuRoot.querySelectorAll('.d4-menu-item')) as HTMLElement[];

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

    // Handle pasting multiple monomers (space / comma / newline separated)
    inputEl.addEventListener('paste', (e: ClipboardEvent) => {
      const pastedText = e.clipboardData?.getData('text');
      if (!pastedText)
        return;

      // Split on newlines first, then parse each line (handles parenthesized symbols like hArg(Et,Et))
      const lines = pastedText.split(/[\n\r]+/).map((l) => l.trim()).filter((l) => l.length > 0);
      const candidates: string[] = [];
      for (const line of lines) {
        const parsed = parseMonomerSymbolList(line);
        // If parseMonomerSymbolList returned a single token but the line has spaces and no commas,
        // split on spaces as well (e.g. "K P F" -> ["K", "P", "F"])
        if (parsed.length === 1 && line.includes(' ') && !line.includes(','))
          candidates.push(...line.split(/\s+/).map((s) => s.trim()).filter((s) => s.length > 0));
        else
          candidates.push(...parsed);
      }

      if (candidates.length <= 1)
        return; // Single symbol: let default paste + autocomplete handle it

      e.preventDefault();

      for (const candidate of candidates) {
        if (allSymbols.includes(candidate) && !selectedMonomers.includes(candidate))
          selectedMonomers.push(candidate);
      }

      renderTags();
      inputEl.value = '';
      hideMenu();
      inputEl.focus();
    });

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

    // --- Bulk add: from monomer library ---
    const librarySection = ui.div([], {style: {marginTop: '8px'}});
    if (libHelper) {
      const libraryInput = ui.input.choice('Add from library', {items: [] as string[], nullable: true}) as DG.ChoiceInput<string>;
      const libraryStatus = ui.divText('', {style: {fontSize: '11px', color: 'var(--green-2)', marginLeft: '8px', minHeight: '16px'}});

      // Load library names asynchronously
      libHelper.getAvaliableLibraryNames().then((libNames) => {
        libraryInput.items = libNames;
      });

      libraryInput.onChanged.subscribe(async () => {
        const libName = libraryInput.value;
        if (!libName)
          return;

        libraryStatus.textContent = 'Loading...';
        try {
          const lib = await libHelper.readSingleLibraryByName(libName);
          if (lib) {
            const symbols = lib.getMonomerSymbolsByType(polymerType);
            const before = selectedMonomers.length;
            addMonomers(symbols);
            const added = selectedMonomers.length - before;
            libraryStatus.textContent = `Added ${added} monomer(s) from "${libName}"`;
          } else {
            libraryStatus.textContent = `Library "${libName}" not found`;
          }
        } catch (err: any) {
          libraryStatus.textContent = `Error loading library`;
        }
        // Clear the selector after use
        libraryInput.value = null as any;
        setTimeout(() => { libraryStatus.textContent = ''; }, 4000);
      });

      librarySection.appendChild(ui.divV([libraryInput.root, libraryStatus]));
    }

    // --- Bulk add: from monomer collection ---
    const collectionSection = ui.div([], {style: {marginTop: '4px'}});
    if (libHelper) {
      const collectionInput = ui.input.choice('Add from collection', {items: [] as string[], nullable: true}) as DG.ChoiceInput<string>;
      const collectionStatus = ui.divText('', {style: {fontSize: '11px', color: 'var(--green-2)', marginLeft: '8px', minHeight: '16px'}});

      // Load collection names asynchronously
      libHelper.listMonomerCollections().then((collectionNames) => {
        collectionInput.items = collectionNames;
      });

      collectionInput.onChanged.subscribe(async () => {
        const collectionName = collectionInput.value;
        if (!collectionName)
          return;

        collectionStatus.textContent = 'Loading...';
        try {
          const collection = await libHelper.readMonomerCollection(collectionName);
          const symbols = collection.monomerSymbols ?? [];
          const before = selectedMonomers.length;
          addMonomers(symbols);
          const added = selectedMonomers.length - before;
          const displayName = collectionName.replace(/\.json$/, '');
          collectionStatus.textContent = `Added ${added} monomer(s) from "${displayName}"`;
        } catch (err: any) {
          collectionStatus.textContent = `Error loading collection`;
        }
        // Clear the selector after use
        collectionInput.value = null as any;
        setTimeout(() => { collectionStatus.textContent = ''; }, 4000);
      });

      collectionSection.appendChild(ui.divV([collectionInput.root, collectionStatus]));
    }

    // --- Clear all button ---
    const clearBtn = ui.button('Clear all', () => {
      selectedMonomers.length = 0;
      renderTags();
    });
    clearBtn.style.fontSize = '11px';
    clearBtn.style.marginTop = '4px';

    renderTags();

    // --- Dialog layout ---
    const contentDiv = ui.div([
      input.root,
      librarySection,
      collectionSection,
      ui.divH([countLabel, clearBtn], {style: {justifyContent: 'space-between', alignItems: 'center'}}),
      tagsHost,
    ], {style: {minWidth: '400px', minHeight: '250px'}});

    const _dlg = ui.dialog({title: 'Select Monomers', showFooter: true})
      .add(contentDiv)
      .onOK(() => { resolve(selectedMonomers); })
      .onCancel(() => { resolve(null); })
      .show({resizable: true});

    inputEl.focus();
    setTimeout(() => { showSuggestions(); }, 0);
  });
}
