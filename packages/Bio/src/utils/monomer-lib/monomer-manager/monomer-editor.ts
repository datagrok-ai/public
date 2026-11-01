/* The editor docked beside the Manage Monomers grid: the library the grid shows; the monomer's structure in
   a sketcher taking the height the fields leave; the fields, scrolling on their own — its symbol, name and
   types two to a row, its R groups as chips that open a menu of caps, its Details — and, pinned at the
   bottom, Save, a line saying what it will do or why it cannot yet, and the Details toggle. With that row
   pinned, the fields grow upward: Details opens above its toggle, an R group's definition under the chips. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer, findProviderWithLibraryName} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType, MonomerType} from '@datagrok-libraries/bio/src/helm/types';
import {MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {INewMonomerForm, MonomerFormFields} from '@datagrok-libraries/bio/src/utils/monomer-ui';
import {mostSimilarNaturalAnalog} from '@datagrok-libraries/bio/src/utils/macromolecule/monomers';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {Control, Form, LabelColumn, Section, ChoiceInput, TextInput, NumberInput, ColorInput,
  MapInput, Dialog, signal, computed, effect, button, statusLine, plural, div, divH, span, iconButton, tableFromMap,
  timestamp, ReadonlySignal, Signal, ChoiceItem, component} from '@datagrok-libraries/u2';

import {MonomerLibManager} from '../lib-manager';
import {RGroupsEditor, rgroupsEditor} from './rgroups-editor';
import {CapGroups} from './cap-groups';
import {MONOMER_DF_COLUMN_NAMES, getCorrectedMolBlock, getCorrectedSmiles, monomerFromDfRow}
  from './monomer-manager';

const POLYMER_TYPES: PolymerType[] =
  [PolymerTypes.PEPTIDE, PolymerTypes.RNA, PolymerTypes.CHEM, PolymerTypes.BLOB, PolymerTypes.G];
const MONOMER_TYPES: MonomerType[] = [MonomerTypes.BACKBONE, MonomerTypes.BRANCH, MonomerTypes.TERMINAL];
/** The narrowest a caption and its input get before the fields go one to a row: a type choice shows its
 * longest item, and the symbol and the name keep room for a dozen characters. */
const COLUMN_WIDTH = 176;
/** The editor height from which Details opens by itself, until the user opens or closes it. */
const DETAILS_OPEN_HEIGHT = 900;

export const DRAW_STRUCTURE = 'Draw the monomer\'s structure';

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

export interface MonomerEditor extends INewMonomerForm {
  readonly library: ChoiceInput;
  readonly sketcher: grok.chem.Sketcher;
  readonly polymerType: ChoiceInput;
  readonly monomerType: ChoiceInput;
  readonly symbol: TextInput;
  readonly monomerName: TextInput;
  readonly id: NumberInput;
  readonly naturalAnalog: TextInput;
  readonly rgroups: Control & RGroupsEditor;
  readonly meta: MapInput;
  readonly colors: {[key in keyof MonomerColors]: ColorInput};
  /** Natural analog, ID, colors and properties: open when the editor is tall enough, until the user
   * opens or closes it themselves. */
  readonly details: Section;
  /** Why the monomer cannot be saved yet, or null. */
  readonly problem: ReadonlySignal<string | null>;
  /** The user changed the monomer, or started a new one, since it was loaded or saved — until then the
   * manager loads the grid's current row into the editor, so moving through the grid shows each monomer. */
  readonly edited: Signal<boolean>;
  /** The grid's library was reloaded: refreshes the list of libraries and what Save would do. */
  libraryReloaded(libraries?: string[], name?: string): void;
  setEmptyMonomer(): void;
  /** The monomer as Save would write it. */
  buildMonomer(): Monomer;
  saveMonomer(): Promise<void>;
  removeMonomers(monomers: Monomer[], libName: string, notify?: boolean, confirm?: boolean): Promise<void>;
}

export function monomerEditor(host: MonomerEditorHost, libraries: string[], libraryName: string):
  Control & MonomerEditor {
  return component((self) => {
    const edited = signal(false);
    /** The structure's own problem: none drawn, no R groups, an R group twice. */
    const structureProblem = signal<string | null>(DRAW_STRUCTURE);
    /** Bumped when the grid's library is reloaded, so "Add" and "Save" follow what it holds. */
    const libraryVersion = signal(0);
    /** The library the grid shows — the Monomer Library input follows it. */
    const shownLibrary = signal(libraryName);
    /** The corrected SMILES the editor put into the sketcher itself — its echo is not an edit. */
    let loadedSmiles: string | null = null;
    /** Set while the editor writes into the sketcher: what it reports meanwhile is not the drawing yet. */
    let settingStructure = false;
    /** The meta of the monomer loaded, its values typed as the library has them. */
    let loadedMeta: NonNullable<Monomer['meta']> = {};
    /** The colors the monomer is drawn with, as loaded, and whether the library sets them itself. */
    let shownColors: MonomerColors = {line: '', background: '', text: ''};
    let ownColors = false;

    const root = div([], 'monomer-manager-editor');
    const labelColumn = new LabelColumn();
    // fill: the library's name takes the row, its arrow and + stay in the dock however narrow
    const library = new ChoiceInput({label: 'Library', name: 'library', items: libraryItems(libraries),
      value: libraryName,
      nullable: false, fill: true, onInput: (name) => name && void host.openLibrary(name)});
    library.follow(() => shownLibrary.value);
    library.addOptions(iconButton('plus', () => host.createLibrary(),
      {tooltip: 'Create a new monomer library'}));

    const sketcher = new DG.chem.Sketcher();
    sketcher.root.classList.add('monomer-manager-sketcher');
    const changes = sketcher.onChanged.subscribe(() => onStructureChanged());
    self.own(() => changes.unsubscribe());
    self.own(() => sketcher.detach());
    // the drawing re-fits when the dock is resized
    const resized = new ResizeObserver(() => sketcher.resize());
    resized.observe(sketcher.root);
    self.own(() => resized.disconnect());

    const polymerType = new ChoiceInput({label: 'Polymer', name: 'polymerType', items: POLYMER_TYPES,
      value: PolymerTypes.PEPTIDE, nullable: false, tooltipText: 'The polymer type the monomer belongs to'});
    const monomerType = new ChoiceInput({label: 'Type', name: 'monomerType', items: MONOMER_TYPES,
      value: MonomerTypes.BACKBONE, nullable: false, tooltipText: 'Monomer type: backbone, branch or terminal'});
    const symbol = new TextInput({label: 'Symbol', name: 'symbol', required: true});
    const monomerName = new TextInput({label: 'Name', name: 'name', required: true});
    const id = new NumberInput({label: 'ID', name: 'id', mode: 'int', value: 0});
    const naturalAnalog = new TextInput({label: 'Analog', name: 'naturalAnalog',
      tooltipText: 'Natural analog: the natural monomer this one is a variant of'});

    const rgroups = rgroupsEditor({canonicalize,
      isValidSmiles: (smiles) => grok.chem.checkSmiles(smiles), onInput: () => edited.value = true});
    const meta = new MapInput({label: 'Properties', name: 'meta'});
    const colorInput = (label: string) => new ColorInput({label, name: label.toLowerCase(), swatchOnly: true});
    const colors = {line: colorInput('Line'), background: colorInput('Background'), text: colorInput('Text')};
    for (const input of [symbol, monomerName, polymerType, monomerType, id, naturalAnalog, meta,
      ...Object.values(colors)])
      input.on('userInput', () => edited.value = true);

    // the first input with a verdict, as Save's tooltip says it: "Symbol: Value can't be empty"
    type Verdict = {label: string, validity: ReadonlySignal<string | null>};
    const firstOf = (inputs: Verdict[]) => {
      for (const input of inputs) {
        const verdict = input.validity.value;
        if (verdict !== null)
          return `${input.label}: ${verdict}`;
      }
      return null;
    };
    const identityProblem = computed(() => firstOf([symbol, monomerName, polymerType, monomerType]));
    // a structure with no R group or one twice is the R groups' to say; nothing drawn is the whole editor's
    const structureOfRGroups = computed(() => {
      const structure = structureProblem.value;
      return structure === DRAW_STRUCTURE ? null : structure;
    });
    const rgroupsProblem = computed(() => structureOfRGroups.value ?? rgroups.problem.value);
    const detailsProblem = computed(() => firstOf([id, naturalAnalog, meta, ...Object.values(colors)]));
    const problem = computed(() => structureProblem.value ?? identityProblem.value ?? rgroupsProblem.value ??
      detailsProblem.value);
    self.effect(() => rgroups.setStructureProblem(structureOfRGroups.value));
    const exists = computed(() => {
      libraryVersion.value;
      const polymer = polymerTypeOf(polymerType.value);
      const monomerSymbol = symbol.value;
      // of this polymer type only: getMonomer falls back across types, and Save adds under the one picked
      return !!polymer && !!monomerSymbol &&
        !!host.activeLib()?.getMonomerSymbolsByType(polymer).includes(monomerSymbol);
    });

    /** What folded Details holds, small enough for the row beside Save: the three colors as dots and how many
     * properties; the analog and the ID in its tooltip. Defaults (no analog, ID 0) are not news. */
    const summary = (): HTMLElement => {
      const count = computed(() => Object.keys(meta.value).filter((key) => key.trim() !== '').length);
      const dots = Object.values(colors).map((color) => {
        const dot = span('', 'monomer-manager-color-dot');
        effect(() => dot.style.background = color.value);
        return dot;
      });
      const props = span(computed(() => count.value === 0 ? '' : String(count.value)), 'monomer-manager-summary-count');
      const shown = div([...dots, props], 'monomer-manager-summary');
      effect(() => {
        const analog = naturalAnalog.value;
        const monomerId = id.value ?? 0;
        shown.title = [analog ? `natural analog ${analog}` : '', monomerId ? `ID ${monomerId}` : '',
          count.value ? plural(count.value, 'property', 'properties') : ''].filter((x) => x).join(' · ') ||
          'Natural analog, ID, colors and properties';
      });
      return shown;
    };

    const getFields = (): MonomerFormFields => ({
      molecule: sketcher.getSmiles() ?? '',
      polymerType: polymerTypeOf(polymerType.peek()) ?? PolymerTypes.PEPTIDE,
      monomerType: monomerTypeOf(monomerType.peek()) ?? MonomerTypes.BACKBONE,
      symbol: symbol.peek(),
      name: monomerName.peek(),
      id: id.peek(),
      naturalAnalog: naturalAnalog.peek(),
    });

    /** Puts `smiles` into the sketcher and answers it corrected, as the drawing will report it. */
    const setStructure = (smiles: string): string => {
      let corrected = '';
      try {
        corrected = smiles ? getCorrectedSmiles([], smiles) : '';
      }
      catch {
        corrected = smiles;
      }
      loadedSmiles = corrected;
      settingStructure = true;
      try {
        sketcher.setSmiles(smiles);
        // the text box above the drawing would still hold what was last typed into it
        sketcher.molInput.value = '';
      } finally {
        settingStructure = false;
      }
      return corrected;
    };

    /** Shows the colors the monomer is drawn with: its own where the library sets them, else the ones the
     * grid uses (its natural analog's, the default grey) — and Save writes them only if they change. */
    const setColors = (monomer: Monomer | null): void => {
      const own = colorsOf(monomer?.meta?.colors);
      ownColors = !!(own.line || own.background || own.text);
      const drawn = !monomer || ownColors ? undefined :
        host.activeLib()?.getMonomerColors(polymerTypeToHelmType(monomer.polymerType), monomer.symbol);
      const hex = (c: string | undefined, fallback: string) => DG.Color.toHtml(DG.Color.fromHtml(c ?? fallback));
      const shown: MonomerColors = {
        line: hex(own.line ?? drawn?.linecolor, '#202020'),
        background: hex(own.background ?? drawn?.backgroundcolor, '#A0A0A0'),
        text: hex(own.text ?? drawn?.textcolor, '#202020'),
      };
      shownColors = shown;
      colors.line.reset(shown.line);
      colors.background.reset(shown.background);
      colors.text.reset(shown.text);
    };

    const checkStructure = (smiles: string): string | null => {
      if (!smiles)
        return DRAW_STRUCTURE;
      const numbers = rgroupNumbers(smiles);
      if (numbers.length === 0) {
        return 'Mark the monomer\'s attachment points as R groups (R1, R2, ...) — the R-group tool in the ' +
          'sketcher, or [*:1] in the SMILES';
      }
      const twice = numbers.find((n, i) => numbers.indexOf(n) !== i);
      return twice === undefined ? null : `R${twice} appears twice in the structure — renumber one of them`;
    };

    const suggestNaturalAnalog = async (smiles: string, polymer: string): Promise<void> => {
      try {
        const capped = CapGroups.capped(smiles, rgroups.rgroups.peek(), canonicalize);
        const analog = await mostSimilarNaturalAnalog(capped, polymer);
        if (analog && !self.isDisposed)
          naturalAnalog.value = analog;
      }
      catch (e) {
        console.error(e);
      }
    };

    /** The R groups, the structure's verdict and the natural analog follow a structure the user drew
     * or a caller set — the drawing's echo of what the editor loaded suggests nothing. */
    const applyStructure = (smiles: string): void => {
      const verdict = checkStructure(smiles);
      structureProblem.value = verdict;
      const numbers = rgroupNumbers(smiles);
      // an R label twice leaves the chips as they were, the R groups saying which one to renumber
      if (new Set(numbers).size !== numbers.length)
        return;
      rgroups.setLabels(numbers.map((n) => `R${n}`));
      // the editor's own load echoes back; only a drawn change suggests a natural analog
      const echo = smiles === loadedSmiles;
      if (!echo) {
        loadedSmiles = null;
        edited.value = true;
      }
      if (!echo && smiles && verdict === null && polymerType.peek())
        void suggestNaturalAnalog(smiles, polymerType.peek()!);
    };

    const onStructureChanged = (): void => {
      if (settingStructure)
        return;
      let smiles = sketcher.getSmiles();
      try {
        smiles = smiles ? getCorrectedSmiles([], smiles) : '';
      }
      catch (e) {
        console.error(e);
      }
      applyStructure(smiles);
    };

    const setFields = (fields: Partial<MonomerFormFields>): void => {
      if (fields.molecule !== undefined) {
        const smiles = setStructure(fields.molecule);
        // a structure set from outside is an edit: its R groups are there at once, not after the drawing lands
        loadedSmiles = null;
        applyStructure(smiles);
        // and the drawing's echo of it, when it lands, is not another one
        loadedSmiles = smiles;
      }
      if (fields.polymerType !== undefined)
        polymerType.value = fields.polymerType;
      if (fields.monomerType !== undefined)
        monomerType.value = fields.monomerType;
      if (fields.symbol !== undefined)
        symbol.value = fields.symbol;
      if (fields.name !== undefined)
        monomerName.value = fields.name;
      if (fields.id !== undefined)
        id.value = fields.id;
      if (fields.naturalAnalog !== undefined)
        naturalAnalog.value = fields.naturalAnalog;
      edited.value = true;
    };

    const setEmptyMonomer = (): void => {
      structureProblem.value = checkStructure(setStructure(''));
      // reset, not set: another monomer starts untouched, with nothing left over from the last one's typing
      symbol.reset('');
      monomerName.reset('');
      id.reset(0);
      naturalAnalog.reset('');
      rgroups.setDefinitions([]);
      loadedMeta = {};
      meta.reset({});
      setColors(null);
      edited.value = false;
    };

    const setMonomer = (monomer: Monomer): void => {
      const smiles = setStructure(monomer.smiles);
      polymerType.reset(monomer.polymerType);
      monomerType.reset(monomer.monomerType);
      symbol.reset(monomer.symbol);
      monomerName.reset(monomer.name);
      id.reset(monomer.id ?? 0);
      naturalAnalog.reset(monomer.naturalAnalog ?? '');
      rgroups.setDefinitions(monomer.rgroups ?? []);
      const loaded = monomer.meta ?? {};
      loadedMeta = loaded;
      meta.reset(Object.fromEntries(Object.keys(loaded)
        .filter((k) => k.toLowerCase() !== 'colors').map((k) => [k, metaText(loaded[k])])));
      structureProblem.value = checkStructure(smiles);
      if (!monomer.naturalAnalog && monomer.polymerType)
        void suggestNaturalAnalog(monomer.smiles, monomer.polymerType);
      setColors(monomer);
      edited.value = false;
    };

    const buildMonomer = (): Monomer => {
      const smiles = getCorrectedSmiles([], sketcher.getSmiles());
      const molfile =
        getCorrectedMolBlock(grok.chem.convert(smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock));
      // a value left as it was keeps its type (true, 1, 124.05), as the library has it
      const written: NonNullable<Monomer['meta']> = {};
      const typed = meta.peek();
      for (const key of Object.keys(typed)) {
        const loaded = loadedMeta[key];
        written[key] = loaded !== undefined && loaded !== null && metaText(loaded) === typed[key] ?
          loaded : typed[key];
      }
      const chosen: MonomerColors = {line: colors.line.peek(), background: colors.background.peek(),
        text: colors.text.peek()};
      const shown = shownColors;
      const changed = chosen.line !== shown.line || chosen.background !== shown.background ||
        chosen.text !== shown.text;
      if (ownColors || changed)
        written.colors = {default: chosen};
      const fields = getFields();
      return {
        symbol: fields.symbol,
        name: fields.name,
        molfile,
        smiles,
        polymerType: fields.polymerType,
        monomerType: fields.monomerType,
        naturalAnalog: fields.naturalAnalog || undefined,
        id: fields.id ?? 0,
        rgroups: rgroups.rgroups.peek(),
        author: DG.User.current().friendlyName,
        createDate: new Date().toISOString(),
        meta: written,
      };
    };

    const addMonomerToLib = async (monomer: Monomer, libName: string): Promise<void> => {
      const df = host.monomersDf();
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
      const same = (i: number) => polymerTypes[i] === monomer.polymerType && symbols[i] === monomer.symbol;
      const sameSymbol = polymerTypes.findIndex((_, i) => same(i));
      // the grid holds canonical SMILES, so the same structure is the same string; a monomer's own does not count
      const sameStructure = structures.findIndex((smi, i) => smi === monomer.smiles && !same(i));

      const save = async () => {
        try {
          const provider = await findProviderWithLibraryName(await host.libManager.getProviders(), libName);
          if (!provider) {
            grok.shell.error(`Library ${libName} not found in any provider`);
            return;
          }
          await provider.updateOrAddMonomersInLibrary(libName, [{...monomer, lib: undefined, wem: undefined}]);
          // saved: the editor follows the grid again, starting with the row the save makes current
          edited.value = false;
          await host.libManager.loadMonomerLib(true);
          await host.refreshTable(monomer.symbol);
          grok.shell.info(`Monomer ${monomer.symbol} was successfully saved in library ${libName}`);
        }
        catch (e) {
          grok.shell.error('Error saving monomer');
          console.error(e);
        }
      };
      // a symbol the library has, and a structure it has under another symbol: each is asked about, both at once
      const prompts: string[] = [];
      const clashes: Monomer[] = [];
      if (sameSymbol >= 0) {
        clashes.push(await monomerFromDfRow(df.rows.get(sameSymbol)));
        prompts.push(`Monomer with symbol '${monomer.symbol}' already exists in library ${libName}.\n` +
          'Are you sure you want to overwrite it?');
      }
      if (sameStructure >= 0) {
        const twin = await monomerFromDfRow(df.rows.get(sameStructure));
        clashes.push(twin);
        prompts.push(`Monomer with the same structure already exists in library ${libName} with different symbol ` +
          `(${twin.symbol}).\nAre you sure you want to duplicate it?`);
      }
      if (clashes.length === 0)
        return await save();
      Dialog.create('Save Monomer', {name: 'saveMonomerDialog'})
        .add(span(prompts.join('\n\n'), 'monomer-manager-prompt'))
        .add(div(clashes.map((m) => monomerCard(m)), 'monomer-manager-cards'))
        .onOK(save, {text: 'Save'})
        .show();
    };

    const saveMonomer = async (): Promise<void> => {
      const verdict = problem.peek();
      if (verdict !== null) {
        grok.shell.warning(verdict);
        return;
      }
      const smiles = sketcher.getSmiles();
      if (!smiles || !grok.chem.checkSmiles(smiles)) {
        grok.shell.warning('Invalid SMILES');
        return;
      }
      const source = host.activeLib()?.source;
      if (!source) {
        grok.shell.warning('Monomer library source is not specified');
        return;
      }
      await addMonomerToLib(buildMonomer(), source);
    };

    const removeMonomers = async (monomers: Monomer[], libName: string, notify = true, confirm = true):
      Promise<void> => {
      const provider = await findProviderWithLibraryName(await host.libManager.getProviders(), libName);
      if (!provider) {
        grok.shell.error(`Library ${libName} not found in any provider`);
        return;
      }
      const several = monomers.length > 1;
      const remove = async () => {
        await provider.deleteMonomersFromLibrary(libName, monomers);
        await host.libManager.loadMonomerLib(true);
        await host.refreshTable();
        // the monomer the editor held is gone: Save must not write it back
        const held = (m: Monomer) => m.symbol === symbol.peek() && m.polymerType === polymerType.peek();
        if (monomers.some(held))
          setEmptyMonomer();
        if (notify) {
          grok.shell.info(`Monomer${several ? 's' : ''} ${monomers.map((m) => m.symbol).join(', ')} ` +
            `${several ? 'were' : 'was'} successfully removed from ${libName} library`);
        }
      };
      if (!confirm)
        return await remove();
      const symbols = monomers.map((m) => m.symbol).join(', ');
      const prompt = several ? `Are you sure you want to remove monomers ${symbols} from ${libName} library?` :
        `Are you sure you want to remove monomer with symbol ${monomers[0].symbol} from ${libName} library?`;
      Dialog.create(`Remove Monomer${several ? 's' : ''}`, {name: 'removeMonomersDialog'})
        .add(span(prompt, 'monomer-manager-prompt'))
        .add(div(monomers.map((m) => monomerCard(m)), 'monomer-manager-cards'))
        .onOK(remove, {text: 'Remove'})
        .show();
    };

    const head = new Form({layout: 'auto', labelColumn}).add(library);
    const fields = new Form({layout: 'auto', labelColumn})
      .addColumns([symbol, monomerName, polymerType, monomerType], {minWidth: COLUMN_WIDTH});
    rgroups.root.classList.add('monomer-manager-rgroups');
    fields.addElement(rgroups.root);
    const save = button(computed(() => exists.value ? 'Save Monomer' : 'Add Monomer'), () => void saveMonomer(),
      {primary: true, disabled: problem});
    // nothing drawn yet is where a new monomer starts, not an error
    const drawing = computed(() => problem.value === DRAW_STRUCTURE);
    const status = statusLine(computed(() => drawing.value ? null : problem.value),
      computed(() => drawing.value ? DRAW_STRUCTURE : `Saves into ${libraryLabel(library.value ?? '')}`));
    status.classList.add('monomer-manager-status');
    const details = new Section({title: 'Details', hint: summary(), validity: detailsProblem,
      chevron: 'end', expanded: false, persistKey: 'bio.monomer-editor.details'});
    details.root.classList.add('monomer-manager-details');
    details.header.classList.add('monomer-manager-details-toggle');
    details.body.classList.add('monomer-manager-details-body');
    // the form's own button row, pinned under the scrolling fields
    const footer = div([], 'monomer-manager-footer');
    fields.addButtons((row) => {
      row.classList.add('monomer-manager-save');
      row.append(save, status, details.header);
    }, {align: 'start', host: footer});
    fields.addElement(details.root)
      .addColumns([naturalAnalog, id], {minWidth: COLUMN_WIDTH, host: details.body})
      .addGroup('Colors', Object.values(colors), {captionsAfter: true, host: details.body})
      .add(meta, details.body);
    head.root.classList.add('monomer-manager-library');
    fields.root.classList.add('monomer-manager-fields');
    // a reversed column keeps the scroll at the bottom: what opens pushes up, toward the sketcher
    const scroller = div([fields.root], 'monomer-manager-scroller');
    // reversed, the scroll runs from 0 at the bottom to minus what is above
    const above = () => scroller.classList.toggle('monomer-manager-more-above',
      scroller.scrollHeight - scroller.clientHeight + scroller.scrollTop > 1);
    self.listen(scroller, 'scroll', above);
    const fit = new ResizeObserver(() => {
      if (!details.userSet && root.clientHeight > 0)
        details.expanded.value = root.clientHeight >= DETAILS_OPEN_HEIGHT;
      above();
    });
    for (const el of [root, scroller, fields.root])
      fit.observe(el);
    self.own(() => fit.disconnect());
    root.append(div([head.root, sketcher.root, scroller, footer], 'monomer-manager-editor-layout'));
    setColors(null);

    return {root, form: root, library, sketcher, polymerType, monomerType, symbol, monomerName, id, naturalAnalog,
      rgroups, meta, colors, details, problem, edited,
      libraryReloaded: (names?: string[], name?: string): void => {
        if (names !== undefined)
          library.setItems(libraryItems(names));
        if (name !== undefined)
          shownLibrary.value = name;
        libraryVersion.value++;
      },
      getFields, setFields, setEmptyMonomer, setMonomer, buildMonomer, saveMonomer, removeMonomers};
  }, {name: 'monomerEditor'});
}

/** Canonical SMILES through Chem, null for what is not one. */
export function canonicalize(smiles: string): string | null {
  try {
    return grok.chem.checkSmiles(smiles) ?
      grok.chem.convert(smiles, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) : null;
  }
  catch {
    return null;
  }
}

export function polymerTypeOf(value: string | null): PolymerType | undefined {
  return POLYMER_TYPES.find((t) => t === value);
}

export function monomerTypeOf(value: string | null): MonomerType | undefined {
  return MONOMER_TYPES.find((t) => t === value);
}

function libraryLabel(name: string): string {
  return name.replace(/\.json$/i, '');
}

function libraryItems(names: string[]): ChoiceItem[] {
  return names.map((value) => ({value, label: libraryLabel(value)}));
}

/** A meta value as the Properties input shows it: an object as its JSON, anything else as its text. */
function metaText(value: string | number | boolean | object | null | undefined): string {
  return value !== null && typeof value === 'object' ? JSON.stringify(value) : String(value ?? '');
}

/** The colors a library keeps in `meta.colors.default`, as an object or as its JSON. */
function colorsOf(colors: string | {default?: Partial<MonomerColors>} | undefined): Partial<MonomerColors> {
  try {
    const parsed: {default?: Partial<MonomerColors>} | undefined =
      typeof colors === 'string' ? JSON.parse(colors) : colors;
    return parsed?.default ?? {};
  }
  catch (e) {
    console.error(e);
    return {};
  }
}

/** The R labels of a SMILES, as written ([1*], [*:1]) and in order. */
function rgroupNumbers(smiles: string): number[] {
  return Array.from(smiles.matchAll(/\[(\d)\*\]|\[\*:(\d)\]/g)).map((m) => Number.parseInt(m[1] ?? m[2]));
}

/** A monomer's structure beside its name, author and date — what a confirmation shows. */
function monomerCard(monomer: Monomer): HTMLElement {
  let image: HTMLElement;
  try {
    image = grok.chem.drawMolecule(getCorrectedSmiles(monomer.rgroups, monomer.smiles, monomer.molfile), 150, 150);
  }
  catch {
    image = span(monomer.symbol);
  }
  const details: Record<string, string | HTMLElement> = {Symbol: monomer.symbol, Name: monomer.name};
  if (monomer.author)
    details['Author'] = monomer.author;
  if (monomer.createDate)
    details['Created'] = timestamp(monomer.createDate);
  return divH([image, tableFromMap(details)], 'monomer-manager-card');
}
