import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {Atom, HelmType, IWebEditorMonomer, GetMonomerResType} from '@datagrok-libraries/bio/src/helm/types';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {
  HELM_REQUIRED_FIELD as REQ
} from '@datagrok-libraries/bio/src/utils/const';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';

import {AmbiguousWebEditorMonomer, GapWebEditorMonomer, getRS, MissingWebEditorMonomer} from './get-monomer-dummy';
import {LibraryWebEditorMonomer} from './get-monomer-of-library';
import {OrgHelmModule, ScilModule} from '../types';
import {RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from '../constants';

import {_package} from '../package';

declare const org: OrgHelmModule;
declare const scil: ScilModule;

const monomerRe = /[\w()]+/;
//** Do not mess with monomer symbol with parenthesis enclosed in square brackets */
const ambMonomerRe = RegExp(String.raw`\(${monomerRe}(,${monomerRe})+\)`);

type GetMonomerOverridingFunc = (
  a: Atom<HelmType> | HelmType, name: string | undefined, monomerLib: IMonomerLib,
  originalGetMonomer: (a: Atom<HelmType> | HelmType, name: string) => GetMonomerResType) => GetMonomerResType;

export function getMonomerOverrideAndLogAlert(
  monomerLib: IMonomerLib, getMonomerOverriding: GetMonomerOverridingFunc, trigger: () => void, logger: ILogger,
): void {
  const logPrefix = 'Helm: overrideGetMonomerAndAlert()';
  const monomers = org.helm.webeditor.Monomers;
  const getMonomerOriginal = monomers.getMonomer.bind(monomers);
  const alertOriginal = scil.Utils.alert;
  try {
    org.helm.webeditor.Monomers.getMonomer =
      (a: Atom<HelmType> | HelmType, name?: string): GetMonomerResType => {
        return getMonomerOverriding(a, name, monomerLib, getMonomerOriginal);
      };
    // Preventing alert message box for missing monomers with compressed Scilligence.JSDraw2.Lite.js
    scil.Utils.alert = (s: string): void => {
      logger.warning(`${logPrefix}, scil.Utils.alert() s = 's'.`);
    };
    trigger();
  } finally {
    monomers.getMonomer = getMonomerOriginal;
    scil.Utils.alert = alertOriginal;
  }
}

/** Inputs logic */
export function getMonomerHandleArgs(
  a: Atom<HelmType> | HelmType, name?: string
): [/** biotype */ HelmType, /** elem */ string] {
  if (!a)
    throw new Error(`Argument 'a' of type Atom or HelmType is mandatory.`);
  let biotype: HelmType;
  let elem: string;
  if ((a as Atom<HelmType>).T === 'ATOM') {
    biotype = (a as Atom<HelmType>).biotype()!;
    elem = (a as Atom<HelmType>).elem;
  } else {
    biotype = a as HelmType;
    elem = org.helm.webeditor.IO.trimBracket(name!);
  }
  return [biotype, elem];
}


/** Substitutes {@link org.helm.webeditor.Monomers.getMonomer()} */
export function getWebEditorMonomer(
  monomerLib: IMonomerLib,
  a: Atom<HelmType> | HelmType, argName?: string,
): IWebEditorMonomer | null {
  const [biotype, elem] = getMonomerHandleArgs(a, argName);
  const pt = helmTypeToPolymerType(biotype);

  /** Get or create {@link Monomer} object (in case it is missing in monomer library current config) */
  let m: Monomer | null = monomerLib.getMonomer(pt, elem);
  if (m && biotype == HelmTypes.LINKER && m[REQ.RGROUPS].length != 2) {
    // Web Editor expects null
    return null;
  }
  if (m && biotype == HelmTypes.SUGAR && m[REQ.RGROUPS].length != 3) {
    // Web Editor expects null
    return null;
  }
  if (!m /* && biotype != HelmTypes.LINKER*/)
    m = monomerLib.addMissingMonomer(pt, elem);

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
      resWem = m.wem = LibraryWebEditorMonomer.fromMonomer(biotype, m);
  }

  return resWem;
}

/** Fills org.helm.webeditor.Monomers dictionary for WebEditor */
export function rewriteLibraries(monomerLib: IMonomerLib): void {
  org.helm.webeditor.Monomers.clear();
  monomerLib!.getPolymerTypes().forEach((polymerType) => {
    const monomerSymbols = monomerLib!.getMonomerSymbolsByType(polymerType);
    monomerSymbols.forEach((monomerSymbol) => {
      let isBroken = false;
      const monomer: Monomer = monomerLib!.getMonomer(polymerType, monomerSymbol)!;
      const webEditorMonomer: IWebEditorMonomer = {
        id: monomerSymbol,
        m: monomer.molfile,
        n: monomer.name,
        na: monomer.naturalAnalog,
        rs: monomer.rgroups.length,
        type: monomer.polymerType,
        mt: monomer.monomerType,
        at: {},
      };

      if (monomer.rgroups.length > 0) {
        // @ts-ignore
        webEditorMonomer.rs = monomer.rgroups.length;
        const at: { [prop: string]: any } = {};
        monomer.rgroups.forEach((it) => {
          at[it[RGROUP_LABEL]] = it[RGROUP_CAP_GROUP_NAME];
        });
        webEditorMonomer.at = at;
      } else if (monomer[SMILES] != null) {
        // @ts-ignore
        webEditorMonomer.rs = Object.keys(getRS(monomer[SMILES].toString())).length;
        webEditorMonomer.at = getRS(monomer[SMILES].toString());
      } else
        isBroken = true;

      if (!isBroken)
        org.helm.webeditor.Monomers.addOneMonomer(webEditorMonomer);
    });
  });

  // Obsolete
  const grid: DG.Grid = grok.shell.tv?.grid;
  if (grid) grid.invalidate();
}
