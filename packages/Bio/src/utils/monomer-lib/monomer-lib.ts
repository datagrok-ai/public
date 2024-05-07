/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {Observable, Subject} from 'rxjs';

import {
  IMonomerLib, Monomer, MonomerType, PolymerType, RGroup
} from '@datagrok-libraries/bio/src/types';
import {
  HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP
} from '@datagrok-libraries/bio/src/utils/const';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

import '../../../css/cell-renderer.css';

/** Wrapper for monomers obtained from different sources. For managing monomere
 * libraries, use MolfileHandler class instead */
export class MonomerLib implements IMonomerLib {
  private _monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

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

  /** Creates missing {@link Monomer} */
  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer {
    let mSet = this._monomers[polymerType];
    if (!mSet)
      mSet = this._monomers[polymerType] = {};
    const m = mSet[monomerSymbol] = {
      [REQ.SYMBOL]: monomerSymbol,
      [REQ.NAME]: monomerSymbol,
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

  getMonomer(polymerType: PolymerType, argMonomerSymbol: string): Monomer | null {
    // Adjust RNA's 'R' for ribose to 'r' and 'P' for phosphate to 'p' for case-sensitive monomer names.
    // There are uppercase 'R' and 'P' at RNA samples in test data 'helm2.csv' but lowercase in HELMCoreLibrary.json
    let monomerSymbol = argMonomerSymbol;
    if (polymerType == 'RNA' && monomerSymbol == 'R')
      monomerSymbol = 'r';
    if (polymerType == 'RNA' && monomerSymbol == 'P')
      monomerSymbol = 'p';

    if (polymerType in this._monomers! && monomerSymbol in this._monomers![polymerType])
      return this._monomers![polymerType][monomerSymbol];
    else
      return null;
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
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: PolymerType, element?: string): string[] {
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
      const rGroupIndices = findAllIndices(molfileHandler.atomTypes, 'R#');
      criterion &&= true;
      return criterion;
    });
    return monomers.map((monomer) => monomer?.symbol!);
  }

  get onChanged(): Observable<any> {
    return this._onChanged;
  }

  private _updateInt(lib: IMonomerLib): void {
    const typesNew = lib.getPolymerTypes();
    const types = this.getPolymerTypes();

    typesNew.forEach((type) => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerSymbolsByType(type);
      monomers.forEach((monomerSymbol) => {
        this._monomers[type][monomerSymbol] = lib.getMonomer(type, monomerSymbol)!;
      });
    });
  }

  public update(lib: IMonomerLib): void {
    this._updateInt(lib);
    this._onChanged.next();
  }

  public updateLibs(libList: IMonomerLib[], reload: boolean = false): void {
    if (reload) this._monomers = {};
    for (const lib of libList) if (!lib.error) this._updateInt(lib);
    this._onChanged.next();
  }

  public clear(): void {
    this._monomers = {};
    this._onChanged.next();
  }

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
      const label = (s: string) => {
        return ui.label(s /* span ? */, {classes: 'ui-input-label'});
      };
      res.append(ui.div([
        label('Name'),
        ui.divText(monomer.name, {classes: 'ui-input-text'})
      ], {classes: 'ui-input-root'}));

      // Structure
      const chemOptions = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      if (monomer.molfile) {
        res.append(ui.div([
          label('Mol file'),
          grok.chem.svgMol(monomer.molfile, undefined, undefined, chemOptions),
        ], {classes: 'ui-input-root'}));
      } else if (monomer.smiles) {
        res.append(ui.div([
          label('Smiles'),
          grok.chem.svgMol(monomer.smiles, undefined, undefined, chemOptions),
        ], {classes: 'ui-input-root'}));
      } else {
        // Unable to get monomer's structure
        res.append(ui.div([
          label('No structure')
        ], {classes: 'ui-input-root'}));
      }

      res.append(ui.div([
        label('Source'),
        ui.divText(monomer.lib?.source ?? 'unknown', {classes: 'ui-input-text'}),
      ], {classes: 'ui-input-root'}));
    } else {
      res.append(ui.divV([
        ui.divText(`Monomer '${monomerSymbol}' of type '${polymerType}' not found.`),
        ui.divText('Open the Context Panel, then expand Manage Libraries'),
      ]));
    }
    return res;
  }
}
