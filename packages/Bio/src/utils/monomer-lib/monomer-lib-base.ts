import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {IMonomerLibBase, Monomer, RGroup} from '@datagrok-libraries/bio/src/types/index';
import {HelmAtom, HelmType, IWebEditorMonomer, MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {getMonomerHandleArgs} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {HelmTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';
import {GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';

import {AmbiguousWebEditorMonomer, GapWebEditorMonomer, MissingWebEditorMonomer} from './web-editor-monomer-dummy';
import {LibraryWebEditorMonomer} from './web-editor-monomer-of-library';

import {_package} from '../../package';
import {Observable, Subject} from 'rxjs';

const monomerRe = /[\w()]+/;
//** Do not mess with monomer symbol with parenthesis enclosed in square brackets */
const ambMonomerRe = RegExp(String.raw`\(${monomerRe}(,${monomerRe})+\)`);

export type MonomerLibDataType = { [polymerType: string]: { [monomerSymbol: string]: Monomer } };

export class MonomerLibBase implements IMonomerLibBase {
  protected _onChanged = new Subject<any>();

  get onChanged(): Observable<any> { return this._onChanged; }

  constructor(
    protected _monomers: MonomerLibDataType,
  ) {}

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

  /** Substitutes {@link org.helm.webeditor.Monomers.getMonomer()} */
  getWebEditorMonomer(a: HelmAtom | HelmType, argName?: string): IWebEditorMonomer | null {
    const [biotype, elem] = getMonomerHandleArgs(a, argName);
    const pt = helmTypeToPolymerType(biotype);

    /** Get or create {@link Monomer} object (in case it is missing in monomer library current config) */
    let m: Monomer | null = this.getMonomer(pt, elem);
    if (m && biotype == HelmTypes.LINKER && m[REQ.RGROUPS].length != 2) {
      // Web Editor expects null
      return null;
    }
    if (m && biotype == HelmTypes.SUGAR && m[REQ.RGROUPS].length != 3) {
      // Web Editor expects null
      return null;
    }
    if (!m /* && biotype != HelmTypes.LINKER*/)
      m = this.addMissingMonomer(pt, elem);

    /** Get or create {@link org,helm.WebEditorMonomer} */
    let resWem: IWebEditorMonomer | null = m.wem ?? null;
    if (!resWem) {
      if (elem === '*')
        resWem = m.wem = new GapWebEditorMonomer(biotype, elem);
      else if (
        (biotype === 'HELM_NUCLETIDE' && elem === 'N') ||
        (biotype === 'HELM_AA' && elem === 'X') ||
        (biotype === 'HELM_CHEM' && false) || // TODO: Ambiguous monomer for CHEM
        ambMonomerRe.test(elem) // e.g. (A,R,_)
      )
        resWem = m.wem = new AmbiguousWebEditorMonomer(biotype, elem);
      else if (!m.lib)
        resWem = m.wem = new MissingWebEditorMonomer(biotype, elem);

      if (!resWem)
        resWem = m.wem = LibraryWebEditorMonomer.fromMonomer(biotype, m, this);
    }

    return resWem!;
  }

  getRS(smiles: string): { [r: string]: string } {
    const newS = smiles.match(/(?<=\[)[^\][]*(?=])/gm);
    const res: { [name: string]: string } = {};
    let el = '';
    let digit;
    if (!!newS) {
      for (let i = 0; i < newS.length; i++) {
        if (newS[i] != null) {
          if (/\d/.test(newS[i])) {
            digit = newS[i][newS[i].length - 1];
            newS[i] = newS[i].replace(/[0-9]/g, '');
            for (let j = 0; j < newS[i].length; j++) {
              if (newS[i][j] != ':')
                el += newS[i][j];
            }
            res['R' + digit] = el;
            el = '';
          }
        }
      }
    }
    return res;
  }
}
