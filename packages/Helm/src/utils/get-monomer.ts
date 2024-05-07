import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as org from 'org';
import scil from 'scil';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {Atom, HelmType, IMonomerLib, Monomer, WebEditorMonomer} from '@datagrok-libraries/bio/src/types';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {PolymerTypes, HelmTypes} from '@datagrok-libraries/bio/src/utils/const';
import {
  HELM_REQUIRED_FIELD as REQ
} from '@datagrok-libraries/bio/src/utils/const';

import {AmbiguousWebEditorMonomer, GapWebEditorMonomer, LibraryWebEditorMonomer} from './dummy-monomer';

import {_package} from '../package';

const monomerRe = /[\w()]+/;
//** Do not mess with monomer symbol with parenthesis enclosed in square brackets */
const ambMonomerRe = RegExp(String.raw`\(${monomerRe}(,${monomerRe})+\)`);

export type GetMonomerResType = WebEditorMonomer | null;

export type GetMonomerFunc = (a: Atom<HelmType> | HelmType, name: string | undefined) => GetMonomerResType;
type GetMonomerOverridingFunc = (
  a: Atom<HelmType> | HelmType, name: string, monomerLib: IMonomerLib,
  originalGetMonomer: (a: Atom<HelmType> | HelmType, name: string) => GetMonomerResType) => GetMonomerResType;

export function getMonomerOverrideAndLogAlert(
  monomerLib: IMonomerLib, getMonomerOverriding: GetMonomerOverridingFunc, trigger: () => void, logger: ILogger,
): void {
  const logPrefix = 'Helm: overrideGetMonomerAndAlert()';
  const monomers = org.helm.webeditor.Monomers;
  const getMonomerOriginal = monomers.getMonomer.bind(monomers);
  const alertOriginal = scil.Utils.alert;
  try {
    org.helm.webeditor.Monomers.getMonomer = (
      a: Atom<HelmType> | HelmType, name: string
    ): GetMonomerResType => {
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
function getMonomerHandleArgs(
  a: Atom<HelmType> | HelmType, name: string
): [/** biotype */ HelmType, /** elem */ string] {
  let s: string;
  let biotype: HelmType;
  if ((a as Atom<HelmType>).T === 'ATOM') {
    biotype = (a as Atom<HelmType>).biotype();
    s = (a as Atom<HelmType>).elem;
  } else {
    biotype = a as HelmType;
    s = org.helm.webeditor.IO.trimBracket(name);
  }
  return [biotype, s];
}


/** Substitutes {@link org.helm.webeditor.Monomers.getMonomer()} */
export function getWebEditorMonomer(
  monomerLib: IMonomerLib,
  a: Atom<HelmType> | HelmType, argName: string,
): WebEditorMonomer | null {
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
  let resWem: WebEditorMonomer | undefined = m.wem;
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

    if (!resWem)
      resWem = m.wem = LibraryWebEditorMonomer.fromMonomer(biotype, m);
  }

  return resWem;
}
