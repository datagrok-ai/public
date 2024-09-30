/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {Observable, Subject} from 'rxjs';

import {
  HelmType, HelmAtom,
  MonomerSetType, MonomerType, PolymerType, IWebEditorMonomer,
} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLibBase, IMonomerLib, IMonomerSet, Monomer, MonomerLibData, MonomerLibSummaryType, RGroup} from '@datagrok-libraries/bio/src/types';
import {HELM_OPTIONAL_FIELDS as OPT, HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {MonomerLibBase, MonomerLibDataType} from './monomer-lib-base';

import {_package} from '../../package';

import '../../../css/cell-renderer.css';

/** Wrapper for monomers obtained from different sources. For managing monomere
 * libraries, use MolfileHandler class instead */
export class MonomerLib extends MonomerLibBase implements IMonomerLib {
  private _duplicateMonomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer[] } } = {};
  public get duplicateMonomers(): { [polymerType: string]: { [monomerSymbol: string]: Monomer[] } } {
    return this._duplicateMonomers;
  }

  private _duplicatesHandled = true;
  public get duplicatesHandled() { return this._duplicatesHandled; }

  private duplicatesNotified: boolean = false;

  constructor(
    monomers: MonomerLibDataType,
    public readonly source: string | undefined = undefined,
    public readonly error: string | undefined = undefined,
  ) {
    super(monomers);
    for (const [_monomerType, monomersOfType] of Object.entries(this._monomers)) {
      for (const [_monomerSymbol, monomer] of Object.entries(monomersOfType))
        monomer.lib = this;
    }
  }

  toJSON(): Monomer[] {
    const resJSON: Monomer[] = [];
    for (const set of Object.values(this._monomers)) {
      for (const m of Object.values(set))
        resJSON.push({...m, lib: undefined, wem: undefined});
    }
    return resJSON;
  }

  getMonomer(polymerType: PolymerType | null, argMonomerSymbol: string): Monomer | null {
    const logPrefix = `Bio: MonomerLib.getMonomer()`;
    // Adjust RNA's 'R' for ribose to 'r' and 'P' for phosphate to 'p' for case-sensitive monomer names.
    // There are uppercase 'R' and 'P' at RNA samples in test data 'helm2.csv' but lowercase in HELMCoreLibrary.json
    let monomerSymbol = argMonomerSymbol;
    if (polymerType == 'RNA' && monomerSymbol == 'R')
      monomerSymbol = 'r';
    if (polymerType == 'RNA' && monomerSymbol == 'P')
      monomerSymbol = 'p';

    let res: Monomer | null = null;

    if (!polymerType) {
      _package.logger.warning(`${logPrefix} symbol '${argMonomerSymbol}', polymerType not specified.`);
      // Assume any polymer type
      for (const [_polymerType, dict] of Object.entries(this._monomers)) {
        res = dict[monomerSymbol];
        if (res) break;
      }
    } else {
      const dict = this._monomers[polymerType];
      res = dict?.[monomerSymbol] ?? null;
    }
    return res;
  }

  private _monomerSets: { [biotype: string /*HelmType*/]: MonomerSetType } | null = null;

  getMonomerSet(biotype: HelmType): MonomerSetType | null {
    const polymerType: PolymerType = helmTypeToPolymerType(biotype);
    if (!this._monomerSets)
      this._monomerSets = {};
    if (!(biotype in this._monomerSets)) {
      for (const [monomerSymbol, monomer] of Object.entries(this._monomers[polymerType])) {

      }
    }
    return this._monomerSets[biotype];
  }

  getPolymerTypes(): PolymerType[] {
    return Object.keys(this._monomers) as PolymerType[];
  }

  getMonomerMolsByPolymerType(polymerType: PolymerType): { [monomerSymbol: string]: string } {
    const res: { [monomerSymbol: string]: string } = {};

    Object.keys(this._monomers[polymerType] ?? {}).forEach((monomerSymbol) => {
      res[monomerSymbol] = this._monomers[polymerType][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerSymbolsByType(polymerType: PolymerType): string[] {
    return Object.keys(this._monomers[polymerType]);
  }

  /** Get a list of monomers with specified element attached to specified
   * R-group
   * WARNING: RGroup numbering starts from 1, not 0*/
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: PolymerType, _element?: string): string[] {
    const monomerSymbols = this.getMonomerSymbolsByType(polymerType);
    let monomers = monomerSymbols.map((sym) => this.getMonomer(polymerType, sym));
    monomers = monomers.filter((el) => el !== null);
    if (monomers.length === 0)
      return [];

    function findAllIndices<T>(arr: T[], element: T): number[] {
      return arr.map((value, index) => (value === element ? index : -1))
        .filter((index) => index !== -1);
    }

    monomers = monomers.filter((monomer) => {
      if (!monomer?.rgroups)
        return false;
      let criterion = monomer?.rgroups.length >= rGroupNumber;
      const molfileHandler = MolfileHandler.getInstance(monomer.molfile);
      const _rGroupIndices = findAllIndices(molfileHandler.atomTypes, 'R#');
      criterion &&= true;
      return criterion;
    });
    return monomers.map((monomer) => monomer?.symbol!);
  }

  private _updateLibInt(lib: IMonomerLib): void {
    const typesNew = lib.getPolymerTypes();
    const types = this.getPolymerTypes();

    typesNew.forEach((type) => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerSymbolsByType(type);
      monomers.forEach((monomerSymbol) => {
        if (this._monomers[type][monomerSymbol]) {
          this._duplicateMonomers[type] ??= {};
          this._duplicateMonomers[type][monomerSymbol] ??= [this._monomers[type][monomerSymbol]];
          this._duplicateMonomers[type][monomerSymbol].push(lib.getMonomer(type, monomerSymbol)!);
        }
        this._monomers[type][monomerSymbol] = lib.getMonomer(type, monomerSymbol)!;
      });
    });
  }

  public update(lib: IMonomerLib): void {
    this._updateLibInt(lib);
    this._onChanged.next();
  }

  public updateLibs(libList: IMonomerLib[], reload: boolean = false): void {
    if (reload)
      this._monomers = {};
    this._duplicateMonomers = {}; // Reset duplicates
    for (const lib of libList)
      if (!lib.error) this._updateLibInt(lib);
    if (Object.entries(this.duplicateMonomers).length > 0) {
      getUserLibSettings().then((settings) => {
        this.assignDuplicatePreferances(settings);
      });
    } else
      this._duplicatesHandled = true;

    this._onChanged.next();
  }

  /** Checks wether all duplicated monomers have set preferences in user settings. overwrites those which have. */
  assignDuplicatePreferances(userSettings: UserLibSettings): boolean {
    let res = true;
    for (const polymerType in this.duplicateMonomers) {
      for (const monomerSymbol in this.duplicateMonomers[polymerType]) {
        if (!userSettings.duplicateMonomerPreferences?.[polymerType]?.[monomerSymbol])
          res = false;
        else {
          const source = userSettings.duplicateMonomerPreferences[polymerType][monomerSymbol];
          const monomer = this.duplicateMonomers[polymerType][monomerSymbol].find((m) => m.lib?.source === source);
          if (!monomer)
            res = false;
          else
            this._monomers[polymerType][monomerSymbol] = monomer;
        }
      }
    }
    this._duplicatesHandled = res;
    return res;
  }

  public clear(): void {
    this._monomers = {};
    this._onChanged.next();
  }

  getSummaryObj(): MonomerLibSummaryType {
    const res: MonomerLibSummaryType = {};
    const ptList: PolymerType[] = this.getPolymerTypes();
    for (const pt of ptList)
      res[pt] = this.getMonomerSymbolsByType(pt).length;
    return res;
  }

  getSummaryDf(): DG.DataFrame {
    const ptList = this.getPolymerTypes();

    const countList: number[] = new Array<number>(ptList.length);
    for (const [pt, i] of wu.enumerate(ptList))
      countList[i] = this.getMonomerSymbolsByType(pt).length;

    const resDf: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('polymerType', ptList),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'count', countList),
    ]);
    return resDf;
  }

  /** @deprecated Keep for backward compatibility */
  getSummary(): string {
    const monTypeList: PolymerType[] = this.getPolymerTypes();
    const resStr: string = monTypeList.length == 0 ? 'empty' : monTypeList.map((monType) => {
      return `${monType} ${this.getMonomerSymbolsByType(monType).length}`;
    }).join('\n');
    return resStr;
  }

  getTooltip(biotype: HelmType, monomerSymbol: string): HTMLElement {
    const polymerType = helmTypeToPolymerType(biotype);
    const res = ui.div([], {classes: 'ui-form ui-tooltip'});
    const monomer = this.getMonomer(polymerType, monomerSymbol);
    if (monomer) {
      // Symbol & Name
      const symbol = monomer[REQ.SYMBOL];
      const _name = monomer[REQ.NAME];
      res.append(ui.divH([
        ui.div([symbol], {style: {fontWeight: 'bolder', textWrap: 'nowrap', marginRight: '6px'}}),
        ui.div([monomer.name])
      ], {style: {display: 'flex', flexDirection: 'row', justifyContent: 'left'}}));

      // Structure
      const chemOptions = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      let structureEl: HTMLElement;
      if (monomer.molfile)
        structureEl = grok.chem.svgMol(monomer.molfile, undefined, undefined, chemOptions);
      else if (monomer.smiles) {
        structureEl = ui.divV([
          grok.chem.svgMol(monomer.smiles, undefined, undefined, chemOptions),
          ui.divText('from smiles', {style: {fontSize: 'smaller'}}),
        ]);
      } else {
        // Unable to get monomer's structure
        structureEl = ui.divText('No structure', {style: {margin: '6px'}});
      }
      res.append(ui.div(structureEl,
        {style: {display: 'flex', flexDirection: 'row', justifyContent: 'center', margin: '6px'}}));

      // Source
      res.append(ui.divText(monomer.lib?.source ?? 'Missed in libraries'));
    } else {
      res.append(ui.divV([
        ui.divText(`Monomer '${monomerSymbol}' of type '${polymerType}' not found.`),
        ui.divText('Open the Context Panel, then expand Manage Libraries'),
      ]));
    }
    return res;
  }

  override(data: MonomerLibData): IMonomerLibBase {
    return new OverriddenMonomerLib(data, this);
  }
}

class OverriddenMonomerLib extends MonomerLibBase {
  constructor(
    private readonly data: MonomerLibData,
    private readonly base: MonomerLibBase
  ) {
    super(data);
  }

  get onChanged(): Observable<any> { return this.base.onChanged; }

  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer {
    return this.base.addMissingMonomer(polymerType, monomerSymbol);
  }

  getMonomer(polymerType: PolymerType | null, monomerSymbol: string): Monomer | null {
    const dataMonomer = this.data[polymerType as string]?.[monomerSymbol];
    const resMonomer = dataMonomer ?? this.base.getMonomer(polymerType, monomerSymbol);
    return resMonomer;
  }
}
