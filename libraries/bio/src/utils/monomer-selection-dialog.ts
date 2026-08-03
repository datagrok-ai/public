/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType, PolymerType} from '../helm/types';
import {IMonomerLib, IMonomerLibHelper, Monomer} from '../types/monomer-library';
import {polymerTypeToHelmType} from './macromolecule/utils';

export function parseMonomerSymbolList(src: string): string[] {
  // L, L-hArg(Et,Et), "hArg(Et,Et)"
  return src.split(/,(?![^(]*\))/)
    .map((s) => {
      s = s.trim();
      if (s.slice(0, 1) === `"` && s.slice(-1) === `"`) s = s.slice(1, -1);
      if (s.slice(0, 1) === `'` && s.slice(-1) === `'`) s = s.slice(1, -1);
      return s.trim();
    })
    .filter((s) => !!s);
}

const MAX_SUGGESTIONS = 20;

/** Embeddable widget for selecting monomers with autocomplete, tag-based display,
 * and bulk-add from monomer libraries or collections.
 * Can be embedded directly into any container (dialog, panel, etc.). */
export class MonomerSelectionWidget {
  private readonly selectedMonomers: string[];
  private readonly monomerLib: IMonomerLib;
  private polymerType: PolymerType;
  private helmType: HelmType;
  private allSymbols: string[];
  private readonly libHelper?: IMonomerLibHelper;

  private readonly tagsHost: HTMLDivElement;
  private readonly countLabel: HTMLDivElement;
  private readonly inputEl: HTMLInputElement;
  private currentMenu: DG.Menu | null = null;
  private menuItems: HTMLElement[] = [];
  private highlightedIndex = -1;

  /** The root element to embed in a container */
  public readonly root: HTMLDivElement;

  constructor(
    monomerLib: IMonomerLib,
    polymerType: PolymerType,
    presetMonomers?: string[],
    libHelper?: IMonomerLibHelper,
  ) {
    this.monomerLib = monomerLib;
    this.polymerType = polymerType;
    this.helmType = polymerTypeToHelmType(polymerType);
    this.allSymbols = monomerLib.getMonomerSymbolsByType(polymerType);
    this.libHelper = libHelper;
    this.selectedMonomers = presetMonomers ? [...presetMonomers] : [];

    // --- Tags display (selected monomers) ---
    this.tagsHost = ui.div([], {style: {display: 'flex', flexWrap: 'wrap', gap: '4px', marginTop: '4px', maxHeight: '150px', overflowY: 'auto'}}) as HTMLDivElement;
    this.countLabel = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-4)', marginTop: '4px'}}) as HTMLDivElement;

    // --- Autocomplete input ---
    const input = ui.input.string('Search', {value: ''});
    this.inputEl = input.input as HTMLInputElement;
    this.inputEl.setAttribute('autocomplete', 'off');
    this.inputEl.placeholder = 'Type to search monomers...';

    this.setupInputListeners();

    // --- Bulk add sections ---
    const librarySection = this.createLibrarySection();
    const collectionSection = this.createCollectionSection();

    // --- Clear all button ---
    const clearBtn = ui.button('Clear all', () => {
      this.selectedMonomers.length = 0;
      this.renderTags();
    });
    clearBtn.style.fontSize = '11px';
    clearBtn.style.marginTop = '4px';

    this.renderTags();

    this.root = ui.div([
      input.root,
      librarySection,
      collectionSection,
      ui.divH([this.countLabel, clearBtn], {style: {justifyContent: 'space-between', alignItems: 'center'}}),
      this.tagsHost,
    ], {style: {minWidth: '300px'}}) as HTMLDivElement;
  }

  /** Returns the currently selected monomer symbols */
  getSelectedMonomers(): string[] {
    return [...this.selectedMonomers];
  }

  /** Sets new monomer selections, replacing existing ones */
  setSelectedMonomers(symbols: string[]): void {
    this.selectedMonomers.length = 0;
    this.selectedMonomers.push(...symbols);
    this.renderTags();
  }

  /** Updates the polymer type and refreshes the available symbols list */
  setPolymerType(polymerType: PolymerType): void {
    this.polymerType = polymerType;
    this.helmType = polymerTypeToHelmType(polymerType);
    this.allSymbols = this.monomerLib.getMonomerSymbolsByType(polymerType);
  }

  /** Focuses the search input */
  focus(): void {
    this.inputEl.focus();
  }

  private updateCountLabel(): void {
    this.countLabel.textContent = this.selectedMonomers.length > 0 ? `${this.selectedMonomers.length} monomer(s) selected` : '';
  }

  private renderTags(): void {
    this.tagsHost.innerHTML = '';
    for (const symbol of this.selectedMonomers) {
      const removeBtn = ui.iconFA('times', () => {
        const idx = this.selectedMonomers.indexOf(symbol);
        if (idx >= 0) {
          this.selectedMonomers.splice(idx, 1);
          this.renderTags();
          this.updateCountLabel();
        }
      }, 'Remove');
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
        const tooltip = this.monomerLib.getTooltip(this.helmType, symbol);
        ui.tooltip.show(tooltip, tag.getBoundingClientRect().left, tag.getBoundingClientRect().bottom + 16);
      });
      tag.addEventListener('mouseleave', () => { ui.tooltip.hide(); });

      this.tagsHost.appendChild(tag);
    }
    this.updateCountLabel();
  }

  private addMonomer(symbol: string): void {
    if (!this.selectedMonomers.includes(symbol)) {
      this.selectedMonomers.push(symbol);
      this.renderTags();
    }
    this.inputEl.value = '';
    this.hideMenu();
    this.inputEl.focus();
  }

  private addMonomers(symbols: string[]): void {
    let added = false;
    for (const symbol of symbols) {
      if (!this.selectedMonomers.includes(symbol)) {
        this.selectedMonomers.push(symbol);
        added = true;
      }
    }
    if (added)
      this.renderTags();
  }

  private hideMenu(): void {
    if (this.currentMenu) {
      this.currentMenu.hide();
      this.currentMenu = null;
    }
    this.menuItems = [];
    this.highlightedIndex = -1;
  }

  private getSuggestions(query: string): {symbol: string, monomer: Monomer | null}[] {
    const q = query.toLowerCase();
    const results: {symbol: string, monomer: Monomer | null, rank: number}[] = [];

    for (const symbol of this.allSymbols) {
      if (this.selectedMonomers.includes(symbol))
        continue;
      const monomer = this.monomerLib.getMonomer(this.polymerType, symbol);
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

  private showSuggestions(): void {
    this.hideMenu();
    const query = this.inputEl.value.trim();
    if (!query)
      return;

    const suggestions = this.getSuggestions(query);
    if (suggestions.length === 0)
      return;

    this.currentMenu = DG.Menu.popup();
    const maxElement = suggestions.reduce((max, s) => {
      const label = s.monomer?.name ? `${s.symbol} - ${s.monomer.name}` : s.symbol;
      if (max.length < label.length)
        return label;
      return max;
    }, '');
    this.currentMenu.item(maxElement, () => {}); // Dummy item to set menu width

    const causedBy = new MouseEvent('mousemove', {clientX: this.inputEl.getBoundingClientRect().left, clientY: this.inputEl.getBoundingClientRect().bottom});
    this.currentMenu.show({causedBy: causedBy,
      y: this.inputEl.offsetHeight + this.inputEl.offsetTop, x: this.inputEl.offsetLeft});

    // collect menu items for keyboard navigation
    setTimeout(() => {
      this.currentMenu?.clear();
      for (const s of suggestions) {
        const label = s.monomer?.name ? `${s.symbol} - ${s.monomer.name}` : s.symbol;
        this.currentMenu?.item(label, () => { this.addMonomer(s.symbol); });
      }

      const menuRoot = document.querySelector('.d4-menu-popup:last-of-type');
      if (menuRoot)
        this.menuItems = Array.from(menuRoot.querySelectorAll('.d4-menu-item')) as HTMLElement[];

      this.highlightedIndex = -1;
    }, 0);
  }

  private updateHighlight(newIndex: number): void {
    if (this.menuItems.length === 0)
      return;
    if (this.highlightedIndex >= 0 && this.highlightedIndex < this.menuItems.length)
      this.menuItems[this.highlightedIndex].classList.remove('d4-menu-item-hover');
    this.highlightedIndex = newIndex;
    if (this.highlightedIndex >= 0 && this.highlightedIndex < this.menuItems.length) {
      this.menuItems[this.highlightedIndex].classList.add('d4-menu-item-hover');
      this.menuItems[this.highlightedIndex].scrollIntoView({block: 'nearest'});
    }
  }

  private setupInputListeners(): void {
    this.inputEl.addEventListener('input', () => { this.showSuggestions(); });

    // Handle pasting multiple monomers (space / comma / newline separated)
    this.inputEl.addEventListener('paste', (e: ClipboardEvent) => {
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
        if (!this.selectedMonomers.includes(candidate))
          this.selectedMonomers.push(candidate);
      }

      this.renderTags();
      this.inputEl.value = '';
      this.hideMenu();
      this.inputEl.focus();
    });

    this.inputEl.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        if (this.menuItems.length > 0) {
          const next = (this.highlightedIndex + 1) % this.menuItems.length;
          this.updateHighlight(next);
        }
      } else if (e.key === 'ArrowUp') {
        e.preventDefault();
        if (this.menuItems.length > 0) {
          const prev = (this.highlightedIndex - 1 + this.menuItems.length) % this.menuItems.length;
          this.updateHighlight(prev);
        }
      } else if (e.key === 'Enter') {
        e.preventDefault();
        e.stopPropagation();
        if (this.highlightedIndex >= 0 && this.highlightedIndex < this.menuItems.length) {
          this.menuItems[this.highlightedIndex].click();
        } else {
          // If input exactly matches a symbol, add it directly
          const val = this.inputEl.value.trim();
          if (val)
            this.addMonomer(val);
        }
      } else if (e.key === 'Escape') {
        this.hideMenu();
      }
    });
  }

  private createLibrarySection(): HTMLDivElement {
    const librarySection = ui.div([], {style: {marginTop: '8px'}}) as HTMLDivElement;
    if (!this.libHelper)
      return librarySection;

    const libHelper = this.libHelper;
    const libraryInput = ui.input.choice('Add from library', {items: [] as string[], nullable: true}) as DG.ChoiceInput<string>;
    const libraryStatus = ui.divText('', {style: {fontSize: '11px', color: 'var(--green-2)', marginLeft: '8px', minHeight: '16px'}});

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
          const symbols = lib.getMonomerSymbolsByType(this.polymerType);
          const before = this.selectedMonomers.length;
          this.addMonomers(symbols);
          const added = this.selectedMonomers.length - before;
          libraryStatus.textContent = `Added ${added} monomer(s) from "${libName}"`;
        } else {
          libraryStatus.textContent = `Library "${libName}" not found`;
        }
      } catch (err: any) {
        libraryStatus.textContent = `Error loading library`;
      }
      libraryInput.value = null as any;
      setTimeout(() => { libraryStatus.textContent = ''; }, 4000);
    });

    librarySection.appendChild(ui.divV([libraryInput.root, libraryStatus]));
    return librarySection;
  }

  private createCollectionSection(): HTMLDivElement {
    const collectionSection = ui.div([], {style: {marginTop: '4px'}}) as HTMLDivElement;
    if (!this.libHelper)
      return collectionSection;

    const libHelper = this.libHelper;
    const collectionInput = ui.input.choice('Add from collection', {items: [] as string[], nullable: true}) as DG.ChoiceInput<string>;
    const collectionStatus = ui.divText('', {style: {fontSize: '11px', color: 'var(--green-2)', marginLeft: '8px', minHeight: '16px'}});

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
        const before = this.selectedMonomers.length;
        this.addMonomers(symbols);
        const added = this.selectedMonomers.length - before;
        const displayName = collectionName.replace(/\.json$/, '');
        collectionStatus.textContent = `Added ${added} monomer(s) from "${displayName}"`;
      } catch (err: any) {
        collectionStatus.textContent = `Error loading collection`;
      }
      collectionInput.value = null as any;
      setTimeout(() => { collectionStatus.textContent = ''; }, 4000);
    });

    collectionSection.appendChild(ui.divV([collectionInput.root, collectionStatus]));
    return collectionSection;
  }
}

/** Shows a dialog for selecting monomers with autocomplete, tag-based display,
 * and bulk-add from monomer libraries or collections.
 * @returns monomer symbols array, or null if cancelled */
export async function showMonomerSelectionDialog(
  monomerLib: IMonomerLib, polymerType: PolymerType, presetMonomers?: string[],
  libHelper?: IMonomerLibHelper,
): Promise<string[] | null> {
  return new Promise<string[] | null>((resolve) => {
    const widget = new MonomerSelectionWidget(monomerLib, polymerType, presetMonomers, libHelper);

    const _dlg = ui.dialog({title: 'Select Monomers', showFooter: true})
      .add(ui.div([widget.root], {style: {minWidth: '400px', minHeight: '250px', maxWidth: '800px'}}))
      .onOK(() => { resolve(widget.getSelectedMonomers()); })
      .onCancel(() => { resolve(null); })
      .show({resizable: true});

    widget.focus();
  });
}
