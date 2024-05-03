/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import * as org from 'org';
import PolymerType = org.helm.PolymerType;
import WebEditorMonomer = org.helm.WebEditorMonomer;

import {Observable, Subject} from 'rxjs';

import {IMonomerLib, Monomer, RGroup} from '@datagrok-libraries/bio/src/types/index';
import {
  HELM_MONOMER_TYPE, HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP
} from '@datagrok-libraries/bio/src/utils/const';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {getRS, helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';

import '../../../css/cell-renderer.css';
import {AmbiguousWebEditorMonomer, DummyWebEditorMonomer, GapWebEditorMonomer} from './dummy-monomer';

const monomerRe = /[\w()]+/;
const ambMonomerRe = RegExp(String.raw`\(${monomerRe}(,${monomerRe})+\)`);

/** Inputs logic */
function getMonomerHandleArgs(
  a: org.helm.IAtom | org.helm.HelmType, name: string
): [/** biotype */ org.helm.HelmType, /** elem */ string] {
  let s: string;
  let biotype: org.helm.HelmType;
  if ((a as org.helm.IAtom).T === 'ATOM') {
    biotype = (a as org.helm.IAtom).biotype();
    s = (a as org.helm.IAtom).elem;
  } else {
    biotype = a as org.helm.HelmType;
    s = org.helm.webeditor.IO.trimBracket(name);
  }
  return [biotype, s];
}

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

  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer {
    const m = this._monomers[polymerType][monomerSymbol] = {
      [REQ.SYMBOL]: monomerSymbol,
      [REQ.NAME]: monomerSymbol,
      [REQ.MOLFILE]: '',
      [REQ.AUTHOR]: 'MISSING',
      [REQ.ID]: -1,
      [REQ.RGROUPS]: [
        /* eslint-disable no-multi-spaces */
        {
          // Samples                     //  PEPTIDE     RNA
          [RGP.CAP_GROUP_SMILES]: '',    // '[*:1][H]'  '[*:1][H]'
          [RGP.ALTERNATE_ID]: '',        // 'R1-H'      'R1-H'
          [RGP.CAP_GROUP_NAME]: '',      // 'H'         'H'
          [RGP.LABEL]: 'R1',             // 'R1'        'R1'
        } as RGroup,
        {
          [RGP.CAP_GROUP_SMILES]: '',    // 'O[*:2]'    '[*:2][H]'
          [RGP.ALTERNATE_ID]: '',        // 'R2-OH'     'R2-H'
          [RGP.CAP_GROUP_NAME]: '',      // 'OH'        'H'
          [RGP.LABEL]: 'R2',             // 'R2'        'R2'
        } as RGroup,
        /* eslint-enable no-multi-spaces */
      ],
      [REQ.SMILES]: '',
      [REQ.POLYMER_TYPE]: polymerType,
      [REQ.MONOMER_TYPE]: undefined as unknown as org.helm.MonomerType, // TODO: Can we get monomerType from atom of POM
      [REQ.CREATE_DATE]: null,
    } as Monomer;
    return m;
  }

  getMonomer(polymerType: string, monomerSymbol: string): Monomer | null {
    if (polymerType in this._monomers! && monomerSymbol in this._monomers![polymerType])
      return this._monomers![polymerType][monomerSymbol];
    else
      return null;
  }

  /** Substitutes {@link org.helm.webeditor.Monomers.getMonomer()} */
  getWebEditorMonomer(
    a: org.helm.IAtom | org.helm.HelmType, name: string
  ): WebEditorMonomer {
    const [biotype, elem] = getMonomerHandleArgs(a, name);
    const pt = helmTypeToPolymerType(biotype);
    let m: Monomer | null = this.getMonomer(pt, name);
    if (!m)
      m = this.addMissingMonomer(pt, name);

    let resWem: WebEditorMonomer | undefined = m.wem;
    if (!resWem) {
      if (name === '*')
        resWem = m.wem = new GapWebEditorMonomer(biotype, elem);
      else if (
        (biotype === 'HELM_NUCLETIDE' && name === 'N') ||
        (biotype === 'HELM_AA' && name === 'X') ||
        (biotype === 'HELM_CHEM' && false) || // TODO: Ambiguous monomer for CHEM
        ambMonomerRe.test(elem)
      )
        resWem = m.wem = new AmbiguousWebEditorMonomer(biotype, elem);

      if (!resWem) {
        resWem = m.wem = new class extends DummyWebEditorMonomer {
          public backgroundcolor: string | undefined = undefined;
          public linecolor: string | undefined = undefined;
          public textcolor: string | undefined = undefined;

          public override readonly at: { [rg: string]: string };

          constructor(biotype: string, dgM: Monomer) {
            super(biotype, dgM.symbol, dgM.name, dgM.molfile);

            const smiles = dgM[REQ.SMILES];
            if (dgM.rgroups.length > 0) {
              const at: { [prop: string]: any } = {};
              dgM.rgroups.forEach((it) => {
                at[it[RGP.LABEL]] = it[RGP.CAP_GROUP_NAME];
              });
              this.at = at;
            } else if (smiles) {
              // Generate R-Groups from SMILES
              this.at = getRS(smiles);
            } else
              throw new Error(`Monomer is broken '${dgM.symbol}' of PolymerType '${dgM[REQ.POLYMER_TYPE]}'.`);
          }
        }(biotype, m);
      }
    }

    return resWem;
  }

  getPolymerTypes(): PolymerType[] {
    return Object.keys(this._monomers) as PolymerType[];
  }

  getMonomerMolsByPolymerType(polymerType: string): { [monomerSymbol: string]: string } {
    const res: { [monomerSymbol: string]: string } = {};

    Object.keys(this._monomers[polymerType] ?? {}).forEach((monomerSymbol) => {
      res[monomerSymbol] = this._monomers[polymerType][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerSymbolsByType(polymerType: string): string[] {
    return Object.keys(this._monomers[polymerType]);
  }

  /** Get a list of monomers with specified element attached to specified
   * R-group
   * WARNING: RGroup numbering starts from 1, not 0*/
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: string, element?: string): string[] {
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
    const monTypeList: string[] = this.getPolymerTypes();
    const resStr: string = monTypeList.length == 0 ? 'empty' : monTypeList.map((monType) => {
      return `${monType} ${this.getMonomerSymbolsByType(monType).length}`;
    }).join('\n');
    return resStr;
  }

  getTooltip(polymerType: string, monomerSymbol: string): HTMLElement {
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
