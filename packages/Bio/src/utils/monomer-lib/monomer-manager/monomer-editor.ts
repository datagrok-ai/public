/* The editor docked beside the Manage Monomers grid: the library the grid shows, the monomer's
   structure in a sketcher, its fields on four tabs and Save, which writes it into that library. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer, findProviderWithLibraryName} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType, MonomerType} from '@datagrok-libraries/bio/src/helm/types';
import {MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {INewMonomerForm, MonomerFormFields} from '@datagrok-libraries/bio/src/utils/monomer-ui';
import {mostSimilarNaturalAnalog} from '@datagrok-libraries/bio/src/utils/macromolecule/monomers';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {Control, Form, TabStrip, ChoiceInput, TextInput, NumberInput, ColorInput, MapInput, Dialog, Tooltip,
  signal, computed, button, div, divH, span, iconButton, tableFromMap, timestamp, Input, ReadonlySignal,
  Signal} from '@datagrok-libraries/u2';

import {MonomerLibManager} from '../lib-manager';
import {RGroupsEditor} from './rgroups-editor';
import {CapGroups} from './cap-groups';
import {MONOMER_DF_COLUMN_NAMES, getCorrectedMolBlock, getCorrectedSmiles, monomerFromDfRow}
  from './monomer-manager';

const POLYMER_TYPES: PolymerType[] =
  [PolymerTypes.PEPTIDE, PolymerTypes.RNA, PolymerTypes.CHEM, PolymerTypes.BLOB, PolymerTypes.G];
const MONOMER_TYPES: MonomerType[] = [MonomerTypes.BACKBONE, MonomerTypes.BRANCH, MonomerTypes.TERMINAL];

interface MonomerColors {
  line: string;
  background: string;
  text: string;
}

export interface MonomerEditorHost {
  readonly libManager: MonomerLibManager;
  /** The library the grid shows. */
  activeLib(): IMonomerLib | null;
  monomersDf(): DG.DataFrame | undefined;
  /** Reloads the grid, scrolling to the monomer with `symbol`. */
  refreshTable(symbol?: string): Promise<void>;
  /** The Monomer Library input picked another library. */
  openLibrary(name: string): Promise<void>;
  createLibrary(): void;
}

export class MonomerEditor extends Control implements INewMonomerForm {
  readonly library: ChoiceInput;
  readonly sketcher: grok.chem.Sketcher;
  readonly polymerType: ChoiceInput;
  readonly monomerType: ChoiceInput;
  readonly symbol: TextInput;
  readonly monomerName: TextInput;
  readonly id: NumberInput;
  readonly naturalAnalog: TextInput;
  readonly rgroups: RGroupsEditor;
  readonly meta: MapInput;
  readonly colors: {[key in keyof MonomerColors]: ColorInput};
  readonly tabs: TabStrip;
  /** Why the monomer cannot be saved yet, or null. */
  readonly problem: ReadonlySignal<string | null>;

  static readonly DRAW = 'Draw the monomer\'s structure';

  /** The structure's own problem: none drawn, no R groups, an R group twice. */
  private readonly _structureProblem = signal<string | null>(MonomerEditor.DRAW);
  /** Bumped when the grid's library is reloaded, so "Add" and "Save" follow what it holds. */
  private readonly _libraryVersion = signal(0);
  /** The library the grid shows — the Monomer Library input follows it. */
  private readonly _libraryName: Signal<string>;
  /** The corrected SMILES the editor put into the sketcher itself — its echo is not an edit. */
  private _loadedSmiles: string | null = null;
  /** Set while the editor writes into the sketcher: what it reports meanwhile is not the drawing yet. */
  private _settingStructure = false;
  /** The meta of the monomer loaded, its values typed as the library has them. */
  private _loadedMeta: NonNullable<Monomer['meta']> = {};
  /** The colors the monomer is drawn with, as loaded, and whether the library sets them itself. */
  private _shownColors: MonomerColors = {line: '', background: '', text: ''};
  private _ownColors = false;

  constructor(private readonly _host: MonomerEditorHost, libraries: string[], libraryName: string) {
    super();
    this._libraryName = signal(libraryName);
    this.root.classList.add('monomer-manager-editor');
    this.name = 'monomerEditor';
    const built = this.runInScope(() => {
      // fill: the library's name takes the row, its arrow and + stay in the dock however narrow
      const library = new ChoiceInput({label: 'Monomer Library', name: 'library', items: libraries, value: libraryName,
        nullable: false, fill: true, onChanged: (name) => name && void this._host.openLibrary(name)});
      library.follow(() => this._libraryName.value);
      library.addOptions(iconButton('plus', () => this._host.createLibrary(),
        {tooltip: 'Create a new monomer library'}));

      const sketcher = new DG.chem.Sketcher();
      sketcher.root.classList.add('monomer-manager-sketcher');
      const changes = sketcher.onChanged.subscribe(() => this._onStructureChanged());
      this.own(() => changes.unsubscribe());
      this.own(() => sketcher.detach());
      // the drawing re-fits when the dock is resized
      const resized = new ResizeObserver(() => sketcher.resize());
      resized.observe(sketcher.root);
      this.own(() => resized.disconnect());

      const fields = {
        polymerType: new ChoiceInput({label: 'Polymer Type', name: 'polymerType', items: POLYMER_TYPES,
          value: PolymerTypes.PEPTIDE, nullable: false}),
        monomerType: new ChoiceInput({label: 'Monomer Type', name: 'monomerType', items: MONOMER_TYPES,
          value: MonomerTypes.BACKBONE, nullable: false}),
        symbol: new TextInput({label: 'Monomer Symbol', name: 'symbol', required: true}),
        monomerName: new TextInput({label: 'Monomer Name', name: 'name', required: true}),
        id: new NumberInput({label: 'Monomer ID', name: 'id', mode: 'int', value: 0}),
        naturalAnalog: new TextInput({label: 'Natural Analog', name: 'naturalAnalog'}),
      };

      const rgroups = new RGroupsEditor({canonicalize: MonomerEditor.canonicalize,
        isValidSmiles: (smiles) => grok.chem.checkSmiles(smiles)});
      const meta = new MapInput({label: 'Properties', name: 'meta'});
      const colorInput = (label: string) => new ColorInput({label, name: label.toLowerCase()});
      const colors = {line: colorInput('Line'), background: colorInput('Background'), text: colorInput('Text')};

      const tabs = new TabStrip({tabs: [
        {id: 'monomer', label: 'Monomer', content: new Form({layout: 'auto'}).addAll(Object.values(fields)).root},
        {id: 'rgroups', label: 'R-groups', content: rgroups.root},
        {id: 'meta', label: 'Meta', content: new Form({layout: 'auto'}).add(meta).root},
        {id: 'colors', label: 'Colors', content: new Form({layout: 'auto'}).addAll(Object.values(colors)).root},
      ]});
      const save = button('Add Monomer', () => void this.saveMonomer(), {primary: true});
      // auto: in a narrow dock the captions go above their inputs instead of squeezing them
      const head = new Form({layout: 'auto'}).add(library);
      // Save stays put under the scrolling body, starting where the editors do, as ui.buttonsInput did
      const footer = new Form({layout: 'normal'}).addButtons((row) => row.append(save), {align: 'editors'});
      return {library, sketcher, ...fields, rgroups, meta, colors, tabs, save, head, footer};
    });
    this.library = built.library;
    this.sketcher = built.sketcher;
    this.polymerType = built.polymerType;
    this.monomerType = built.monomerType;
    this.symbol = built.symbol;
    this.monomerName = built.monomerName;
    this.id = built.id;
    this.naturalAnalog = built.naturalAnalog;
    this.rgroups = built.rgroups;
    this.meta = built.meta;
    this.colors = built.colors;
    this.tabs = built.tabs;
    const save = built.save;

    // each tab answers for its own inputs, so the tab with the problem is the one marked
    type Verdict = {label: string, validity: ReadonlySignal<string | null>};
    const firstOf = (inputs: Verdict[]) => {
      for (const input of inputs) {
        const verdict = input.validity.value;
        if (verdict !== null)
          return `${input.label}: ${verdict}`;
      }
      return null;
    };
    const monomerRequired = computed(() => {
      for (const input of [this.polymerType, this.monomerType, this.symbol, this.monomerName]) {
        if (Input.isEmpty(input.value.value))
          return `${input.label} is required`;
      }
      return null;
    });
    // what was typed wrong; an empty required field is monomerRequired's to say
    const monomerVerdicts = computed(() => firstOf([this.id, this.naturalAnalog]));
    const monomerTab = computed(() => monomerRequired.value ?? monomerVerdicts.value);
    // a structure with no R group or one twice is the R-groups tab's to say; nothing drawn is no tab's
    const structureOfRGroups = computed(() => {
      const structure = this._structureProblem.value;
      return structure === MonomerEditor.DRAW ? null : structure;
    });
    const rgroupsTab = computed(() => structureOfRGroups.value ?? this.rgroups.problem.value);
    const metaTab = computed(() => firstOf([this.meta]));
    const colorsTab = computed(() => firstOf(Object.values(this.colors)));
    this.problem = computed(() => this._structureProblem.value ?? monomerTab.value ?? rgroupsTab.value ??
      metaTab.value ?? colorsTab.value);
    this.effect(() => {
      // a new monomer's empty symbol and name mark the tab only once there is something to name
      this.tabs.setInvalid('monomer', this._structureProblem.value === MonomerEditor.DRAW ?
        monomerVerdicts.value : monomerTab.value);
      this.tabs.setInvalid('rgroups', rgroupsTab.value);
      this.tabs.setInvalid('meta', metaTab.value);
      this.tabs.setInvalid('colors', colorsTab.value);
      this.rgroups.setStructureProblem(structureOfRGroups.value);
    });
    const exists = computed(() => {
      this._libraryVersion.value;
      const polymerType = MonomerEditor.polymerTypeOf(this.polymerType.value.value);
      const symbol = this.symbol.value.value;
      // of this polymer type only: getMonomer falls back across types, and Save adds under the one picked
      return !!polymerType && !!symbol &&
        !!this._host.activeLib()?.getMonomerSymbolsByType(polymerType).includes(symbol);
    });
    this.effect(() => {
      save.disabled = this.problem.value !== null;
      save.textContent = exists.value ? 'Save Monomer' : 'Add Monomer';
    });
    Tooltip.bind(save, () => this.problem.peek() ?? `Save the monomer into ${this.library.value.peek()}`, this.scope);

    built.head.root.classList.add('monomer-manager-library');
    built.footer.root.classList.add('monomer-manager-save');
    this.root.append(div([built.head.root, this.sketcher.root, this.tabs.root], 'monomer-manager-editor-body'),
      built.footer.root);
    this._setColors(null);
  }

  /** Canonical SMILES through Chem, null for what is not one. */
  static canonicalize(smiles: string): string | null {
    try {
      return grok.chem.checkSmiles(smiles) ?
        grok.chem.convert(smiles, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) : null;
    } catch {
      return null;
    }
  }

  static polymerTypeOf(value: string | null): PolymerType | undefined {
    return POLYMER_TYPES.find((t) => t === value);
  }

  static monomerTypeOf(value: string | null): MonomerType | undefined {
    return MONOMER_TYPES.find((t) => t === value);
  }

  get form(): HTMLElement {
    return this.root;
  }

  /** The grid's library was reloaded: refreshes the list of libraries and what Save would do. */
  libraryReloaded(libraries?: string[], name?: string): void {
    if (libraries !== undefined)
      this.library.setItems(libraries);
    if (name !== undefined)
      this._libraryName.value = name;
    this._libraryVersion.value++;
  }

  getFields(): MonomerFormFields {
    return {
      molecule: this.sketcher.getSmiles() ?? '',
      polymerType: MonomerEditor.polymerTypeOf(this.polymerType.value.peek()) ?? PolymerTypes.PEPTIDE,
      monomerType: MonomerEditor.monomerTypeOf(this.monomerType.value.peek()) ?? MonomerTypes.BACKBONE,
      symbol: this.symbol.value.peek(),
      name: this.monomerName.value.peek(),
      id: this.id.value.peek(),
      naturalAnalog: this.naturalAnalog.value.peek(),
    };
  }

  setFields(fields: Partial<MonomerFormFields>): void {
    if (fields.molecule !== undefined) {
      const smiles = this._setStructure(fields.molecule);
      // a structure set from outside is an edit: its R groups are there at once, not after the drawing lands
      this._loadedSmiles = null;
      this._applyStructure(smiles);
      // and the drawing's echo of it, when it lands, is not another one
      this._loadedSmiles = smiles;
    }
    if (fields.polymerType !== undefined)
      this.polymerType.value.value = fields.polymerType;
    if (fields.monomerType !== undefined)
      this.monomerType.value.value = fields.monomerType;
    if (fields.symbol !== undefined)
      this.symbol.value.value = fields.symbol;
    if (fields.name !== undefined)
      this.monomerName.value.value = fields.name;
    if (fields.id !== undefined)
      this.id.value.value = fields.id;
    if (fields.naturalAnalog !== undefined)
      this.naturalAnalog.value.value = fields.naturalAnalog;
  }

  setEmptyMonomer(): void {
    this._structureProblem.value = this._checkStructure(this._setStructure(''));
    // reset, not set: another monomer starts untouched, with nothing left over from the last one's typing
    this.symbol.reset('');
    this.monomerName.reset('');
    this.id.reset(0);
    this.naturalAnalog.reset('');
    this.rgroups.setDefinitions([]);
    this._loadedMeta = {};
    this.meta.reset({});
    this._setColors(null);
  }

  setMonomer(monomer: Monomer): void {
    const smiles = this._setStructure(monomer.smiles);
    this.polymerType.reset(monomer.polymerType);
    this.monomerType.reset(monomer.monomerType);
    this.symbol.reset(monomer.symbol);
    this.monomerName.reset(monomer.name);
    this.id.reset(monomer.id ?? 0);
    this.naturalAnalog.reset(monomer.naturalAnalog ?? '');
    this.rgroups.setDefinitions(monomer.rgroups ?? []);
    const meta = monomer.meta ?? {};
    this._loadedMeta = meta;
    this.meta.reset(Object.fromEntries(Object.keys(meta)
      .filter((k) => k.toLowerCase() !== 'colors').map((k) => [k, MonomerEditor.metaText(meta[k])])));
    this._structureProblem.value = this._checkStructure(smiles);
    if (!monomer.naturalAnalog && monomer.polymerType)
      void this._suggestNaturalAnalog(monomer.smiles, monomer.polymerType);
    this._setColors(monomer);
  }

  /** A meta value as the Meta tab shows it: an object as its JSON, anything else as its text. */
  static metaText(value: string | number | boolean | object | null | undefined): string {
    return value !== null && typeof value === 'object' ? JSON.stringify(value) : String(value ?? '');
  }

  /** The colors a library keeps in `meta.colors.default`, as an object or as its JSON. */
  static colorsOf(colors: string | {default?: Partial<MonomerColors>} | undefined): Partial<MonomerColors> {
    try {
      const parsed: {default?: Partial<MonomerColors>} | undefined =
        typeof colors === 'string' ? JSON.parse(colors) : colors;
      return parsed?.default ?? {};
    } catch (e) {
      console.error(e);
      return {};
    }
  }

  /** The monomer as Save would write it. */
  buildMonomer(): Monomer {
    const smiles = getCorrectedSmiles([], this.sketcher.getSmiles());
    const molfile = getCorrectedMolBlock(grok.chem.convert(smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock));
    // a value left as it was keeps its type (true, 1, 124.05), as the library has it
    const meta: NonNullable<Monomer['meta']> = {};
    const edited = this.meta.value.peek();
    for (const key of Object.keys(edited)) {
      const loaded = this._loadedMeta[key];
      meta[key] = loaded !== undefined && loaded !== null && MonomerEditor.metaText(loaded) === edited[key] ?
        loaded : edited[key];
    }
    const colors: MonomerColors = {line: this.colors.line.value.peek(), background: this.colors.background.value.peek(),
      text: this.colors.text.value.peek()};
    const shown = this._shownColors;
    const changed = colors.line !== shown.line || colors.background !== shown.background || colors.text !== shown.text;
    if (this._ownColors || changed)
      meta.colors = {default: colors};
    const fields = this.getFields();
    return {
      symbol: fields.symbol,
      name: fields.name,
      molfile,
      smiles,
      polymerType: fields.polymerType,
      monomerType: fields.monomerType,
      naturalAnalog: fields.naturalAnalog || undefined,
      id: fields.id ?? 0,
      rgroups: this.rgroups.rgroups.peek(),
      author: DG.User.current().friendlyName,
      createDate: new Date().toISOString(),
      meta,
    };
  }

  async saveMonomer(): Promise<void> {
    const problem = this.problem.peek();
    if (problem !== null) {
      grok.shell.warning(problem);
      return;
    }
    const smiles = this.sketcher.getSmiles();
    if (!smiles || !grok.chem.checkSmiles(smiles)) {
      grok.shell.warning('Invalid SMILES');
      return;
    }
    const source = this._host.activeLib()?.source;
    if (!source) {
      grok.shell.warning('Monomer library source is not specified');
      return;
    }
    await this._addMonomerToLib(this.buildMonomer(), source);
  }

  async removeMonomers(monomers: Monomer[], libName: string, notify = true, confirm = true): Promise<void> {
    const provider = await findProviderWithLibraryName(await this._host.libManager.getProviders(), libName);
    if (!provider) {
      grok.shell.error(`Library ${libName} not found in any provider`);
      return;
    }
    const plural = monomers.length > 1;
    const remove = async () => {
      await provider.deleteMonomersFromLibrary(libName, monomers);
      await this._host.libManager.loadMonomerLib(true);
      await this._host.refreshTable();
      // the monomer the editor held is gone: Save must not write it back
      const edited = (m: Monomer) =>
        m.symbol === this.symbol.value.peek() && m.polymerType === this.polymerType.value.peek();
      if (monomers.some(edited))
        this.setEmptyMonomer();
      if (notify) {
        grok.shell.info(`Monomer${plural ? 's' : ''} ${monomers.map((m) => m.symbol).join(', ')} ` +
          `${plural ? 'were' : 'was'} successfully removed from ${libName} library`);
      }
    };
    if (!confirm)
      return await remove();
    const prompt = plural ?
      `Are you sure you want to remove monomers ${monomers.map((m) => m.symbol).join(', ')} from ${libName} library?` :
      `Are you sure you want to remove monomer with symbol ${monomers[0].symbol} from ${libName} library?`;
    Dialog.create(`Remove Monomer${plural ? 's' : ''}`, {name: 'removeMonomersDialog'})
      .add(span(prompt, 'monomer-manager-prompt'))
      .add(div(monomers.map((m) => MonomerEditor.monomerCard(m)), 'monomer-manager-cards'))
      .onOK(remove, {text: 'Remove'})
      .show();
  }

  /** A monomer's structure beside its name, author and date — what a confirmation shows. */
  static monomerCard(monomer: Monomer): HTMLElement {
    let image: HTMLElement;
    try {
      image = grok.chem.drawMolecule(getCorrectedSmiles(monomer.rgroups, monomer.smiles, monomer.molfile), 150, 150);
    } catch {
      image = span(monomer.symbol);
    }
    const details: Record<string, string | HTMLElement> = {Symbol: monomer.symbol, Name: monomer.name};
    if (monomer.author)
      details['Author'] = monomer.author;
    if (monomer.createDate)
      details['Created'] = timestamp(monomer.createDate);
    return divH([image, tableFromMap(details)], 'monomer-manager-card');
  }

  /** Puts `smiles` into the sketcher and answers it corrected, as the drawing will report it. */
  private _setStructure(smiles: string): string {
    let corrected = '';
    try {
      corrected = smiles ? getCorrectedSmiles([], smiles) : '';
    } catch {
      corrected = smiles;
    }
    this._loadedSmiles = corrected;
    this._settingStructure = true;
    try {
      this.sketcher.setSmiles(smiles);
      // the text box above the drawing would still hold what was last typed into it
      this.sketcher.molInput.value = '';
    } finally {
      this._settingStructure = false;
    }
    return corrected;
  }

  /** Shows the colors the monomer is drawn with: its own where the library sets them, else the ones the
   * grid uses (its natural analog's, the default grey) — and Save writes them only if they change. */
  private _setColors(monomer: Monomer | null): void {
    const own = MonomerEditor.colorsOf(monomer?.meta?.colors);
    this._ownColors = !!(own.line || own.background || own.text);
    const drawn = !monomer || this._ownColors ? undefined :
      this._host.activeLib()?.getMonomerColors(polymerTypeToHelmType(monomer.polymerType), monomer.symbol);
    const hex = (c: string | undefined, fallback: string) => DG.Color.toHtml(DG.Color.fromHtml(c ?? fallback));
    const shown: MonomerColors = {
      line: hex(own.line ?? drawn?.linecolor, '#202020'),
      background: hex(own.background ?? drawn?.backgroundcolor, '#A0A0A0'),
      text: hex(own.text ?? drawn?.textcolor, '#202020'),
    };
    this._shownColors = shown;
    this.colors.line.reset(shown.line);
    this.colors.background.reset(shown.background);
    this.colors.text.reset(shown.text);
  }

  /** The R labels of a SMILES, as written ([1*], [*:1]) and in order. */
  static rgroupNumbers(smiles: string): number[] {
    return Array.from(smiles.matchAll(/\[(\d)\*\]|\[\*:(\d)\]/g)).map((m) => Number.parseInt(m[1] ?? m[2]));
  }

  private _checkStructure(smiles: string): string | null {
    if (!smiles)
      return MonomerEditor.DRAW;
    const numbers = MonomerEditor.rgroupNumbers(smiles);
    if (numbers.length === 0)
      return 'Mark the monomer\'s attachment points as R groups (R1, R2, ...)';
    const twice = numbers.find((n, i) => numbers.indexOf(n) !== i);
    return twice === undefined ? null : `R${twice} appears twice in the structure — renumber one of them`;
  }

  private _onStructureChanged(): void {
    if (this._settingStructure)
      return;
    let smiles = this.sketcher.getSmiles();
    try {
      smiles = smiles ? getCorrectedSmiles([], smiles) : '';
    } catch (e) {
      console.error(e);
    }
    this._applyStructure(smiles);
  }

  /** The R groups, the structure's verdict and the natural analog follow a structure the user drew
   * or a caller set — the drawing's echo of what the editor loaded suggests nothing. */
  private _applyStructure(smiles: string): void {
    const problem = this._checkStructure(smiles);
    this._structureProblem.value = problem;
    const numbers = MonomerEditor.rgroupNumbers(smiles);
    // an R label twice leaves the rows as they were, the tab saying which one to renumber
    if (new Set(numbers).size !== numbers.length)
      return;
    this.rgroups.setLabels(numbers.map((n) => `R${n}`));
    // the editor's own load echoes back; only a drawn change suggests a natural analog
    const echo = smiles === this._loadedSmiles;
    if (!echo)
      this._loadedSmiles = null;
    if (!echo && smiles && problem === null && this.polymerType.value.peek())
      void this._suggestNaturalAnalog(smiles, this.polymerType.value.peek()!);
  }

  private async _suggestNaturalAnalog(smiles: string, polymerType: string): Promise<void> {
    try {
      const capped = CapGroups.capped(smiles, this.rgroups.rgroups.peek(), MonomerEditor.canonicalize);
      const analog = await mostSimilarNaturalAnalog(capped, polymerType);
      if (analog && !this.scope.isDisposed)
        this.naturalAnalog.value.value = analog;
    } catch (e) {
      console.error(e);
    }
  }

  private async _addMonomerToLib(monomer: Monomer, libName: string): Promise<void> {
    const df = this._host.monomersDf();
    if (!df) {
      grok.shell.error('No monomer library loaded');
      return;
    }
    const polymerTypes = df.col(MONOMER_DF_COLUMN_NAMES.POLYMER_TYPE)?.toList() ?? [];
    const symbols = df.col(MONOMER_DF_COLUMN_NAMES.SYMBOL)?.toList() ?? [];
    const structures = df.col(MONOMER_DF_COLUMN_NAMES.MONOMER)?.toList() ?? [];
    if (polymerTypes.length !== symbols.length || polymerTypes.length !== structures.length) {
      grok.shell.error('Monomer library data frame is corrupted');
      return;
    }
    const sameSymbol = polymerTypes.findIndex((pt, i) => pt === monomer.polymerType && symbols[i] === monomer.symbol);
    // the grid holds canonical SMILES, so the same structure is the same string
    const sameStructure = structures.findIndex((smi) => smi === monomer.smiles);

    const save = async () => {
      try {
        const provider = await findProviderWithLibraryName(await this._host.libManager.getProviders(), libName);
        if (!provider) {
          grok.shell.error(`Library ${libName} not found in any provider`);
          return;
        }
        await provider.updateOrAddMonomersInLibrary(libName, [{...monomer, lib: undefined, wem: undefined}]);
        await this._host.libManager.loadMonomerLib(true);
        await this._host.refreshTable(monomer.symbol);
        grok.shell.info(`Monomer ${monomer.symbol} was successfully saved in library ${libName}`);
      } catch (e) {
        grok.shell.error('Error saving monomer');
        console.error(e);
      }
    };
    let clash: Monomer | null = null;
    let prompt = '';
    if (sameSymbol >= 0) {
      clash = await monomerFromDfRow(df.rows.get(sameSymbol));
      prompt = `Monomer with symbol '${monomer.symbol}' already exists in library ${libName}.\n` +
        'Are you sure you want to overwrite it?';
    } else if (sameStructure >= 0) {
      clash = await monomerFromDfRow(df.rows.get(sameStructure));
      prompt = `Monomer with the same structure already exists in library ${libName} with different symbol ` +
        `(${clash.symbol}).\nAre you sure you want to duplicate it?`;
    }
    if (clash === null)
      return await save();
    Dialog.create('Save Monomer', {name: 'saveMonomerDialog'})
      .add(span(prompt, 'monomer-manager-prompt'))
      .add(MonomerEditor.monomerCard(clash))
      .onOK(save, {text: 'Save'})
      .show();
  }
}
