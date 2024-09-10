/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {Observable, Subject} from 'rxjs';

import {HelmType, MonomerSetType, MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerSet, Monomer, MonomerLibSummaryType, RGroup} from '@datagrok-libraries/bio/src/types';
import {HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {GapOriginals} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';

import '../../../css/cell-renderer.css';

import {_package} from '../../package';
import {getUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

/** Wrapper for monomers obtained from different sources. For managing monomere
 * libraries, use MolfileHandler class instead */
export class MonomerLib implements IMonomerLib {
  private _monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();
  private _duplicateMonomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer[] } } = {};
  public get duplicateMonomers(): { [polymerType: string]: { [monomerSymbol: string]: Monomer[] } } {
    return this._duplicateMonomers;
  }
  private _duplicatesHandled = true;
  public get duplicatesHandled() { return this._duplicatesHandled; }
  private duplicatesNotified: boolean = false;

  constructor(
    monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } },
    public readonly source: string | undefined = undefined,
    public readonly error: string | undefined = undefined,
  ) {
    this._monomers = monomers;
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

  /** Creates missing {@link Monomer} */
  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer {
    let mSet = this._monomers[polymerType];
    if (!mSet)
      mSet = this._monomers[polymerType] = {};

    let monomerName: string = monomerSymbol;
    if (monomerSymbol === GapOriginals[NOTATION.HELM])
      monomerName = 'Gap';
    else if (polymerType === PolymerTypes.PEPTIDE && monomerSymbol === 'X')
      monomerName = 'Any';
    else if (polymerType === PolymerTypes.RNA && monomerSymbol === 'N')
      monomerName = 'Any';

    const m = mSet[monomerSymbol] = {
      [REQ.SYMBOL]: monomerSymbol,
      [REQ.NAME]: monomerName,
      [REQ.MOLFILE]: '',
      [REQ.AUTHOR]: 'MISSING',
      [REQ.ID]: -1,
      [REQ.RGROUPS]:
        wu.count(1).take(9).map((i) => {
          return {
            /* eslint-disable no-multi-spaces */
            // Samples                        //  PEPTIDE     RNA
            [RGP.CAP_GROUP_SMILES]: '',       // '[*:1][H]'  '[*:1][H]'
            [RGP.ALTERNATE_ID]: '',           // 'R1-H'      'R1-H'
            [RGP.CAP_GROUP_NAME]: '',         // 'H'         'H'
            [RGP.LABEL]: `R${i.toString()}`,  // 'R1'        'R1'
            /* eslint-enable no-multi-spaces */
          } as RGroup;
        }).toArray(),
      [REQ.SMILES]: '',
      [REQ.POLYMER_TYPE]: polymerType,
      [REQ.MONOMER_TYPE]: undefined as unknown as MonomerType, // TODO: Can we get monomerType from atom of POM
      [REQ.CREATE_DATE]: null,
    } as Monomer;
    return m;
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
      res = dict ? dict[monomerSymbol] : null;
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

  get onChanged(): Observable<any> {
    return this._onChanged;
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

  getTooltip(polymerType: PolymerType, monomerSymbol: string): HTMLElement {
    // getTooltip(monomer: Monomer): HTMLElement;
    // getTooltip(monomerOrPolymerType: string | Monomer, symbol?: string): HTMLElement {
    //   let polymerType: string;
    //   let monomerSymbol: string;
    //   if (typeof monomerOrPolymerType === 'string' || monomerOrPolymerType instanceof String) {
    //     polymerType = monomerOrPolymerType as string;
    //     monomerSymbol = symbol!;
    //   } else {
    //     const m = monomerOrPolymerType as Monomer;
    //     polymerType = m[HELM_REQUIRED_FIELD.POLYMER_TYPE];
    //     monomerSymbol = m[HELM_REQUIRED_FIELD.SYMBOL];
    //   }
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
}
