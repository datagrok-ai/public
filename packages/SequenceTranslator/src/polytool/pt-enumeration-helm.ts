import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {HelmType, JSDraw2ModuleType, OrgType, HelmAtom, HelmMol} from '@datagrok-libraries/bio/src/helm/types';

import {Chain} from './pt-conversion';
import {getAvailableMonomers} from './utils';
import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes, PolyToolPlaceholders} from './types';

export const PT_HELM_EXAMPLE = 'PEPTIDE1{[R].[F].[T].[G].[H].[F].[G].[A].[A].[Y].[P].[E].[NH2]}$$$$';

/** Initialized by getHelmHelper via init Helm package */
declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;

function polyToolEnumeratorCore(m: HelmMol, position: number, monomerList: string[]): HelmMol[] {
  const resMolList: HelmMol[] = new Array<HelmMol>(monomerList.length);
  for (let i = 0; i < monomerList.length; i++) {
    const symbolName = monomerList[i];
    const resM = resMolList[i] = m.clone();
    resM.atoms[position].elem = symbolName;
  }
  return resMolList;
}

/**
 * @param {string} helm Molecule string Helm format
 * @param  placeholders Placeholders by zero-based position key
 * @returns {string[]} List of enumerated molecules in Helm format
 */
function getPtEnumeratorSingle(helm: string, placeholders: PolyToolPlaceholders): string[] {
  const molHandler = new JSDraw2.MolHandler<HelmType>();
  const plugin = new org.helm.webeditor.Plugin(molHandler);
  const io = org.helm.webeditor.IO;

  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, helm, origin, undefined);

  const coreResList: HelmMol[][] = Object.entries(placeholders)
    .map(([p, monomerList]: [string, string[]]) => polyToolEnumeratorCore(molHandler.m, parseInt(p), monomerList));
  const resMolList = coreResList.reduce((acc, posList) => acc.concat(posList), []);

  const resHelmList = resMolList.map((m: HelmMol) => org.helm.webeditor.IO.getHelm(m)!);
  return resHelmList;
}

function getPtEnumeratorMatrix(helm: string, placeholders: PolyToolPlaceholders): string[] {
  const molHandler = new JSDraw2.MolHandler<HelmType>();
  const plugin = new org.helm.webeditor.Plugin(molHandler);
  const io = org.helm.webeditor.IO;

  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, helm, origin, undefined);

  let resMolList = [molHandler.m];
  for (const [p, monomerList] of Object.entries(placeholders)) {
    const pos: number = parseInt(p);
    const posResMolList: HelmMol[][] = resMolList.map((m: HelmMol) => polyToolEnumeratorCore(m, pos, monomerList));
    resMolList = posResMolList.reduce((acc, l) => acc.concat(l), []);
  }

  const resHelmList = resMolList.map((m: HelmMol) => org.helm.webeditor.IO.getHelm(m)!);
  return resHelmList;
}

export function getPtEnumeratorHelm(helm: string, params: PolyToolEnumeratorParams): string[] {
  switch (params.type) {
  case PolyToolEnumeratorTypes.Single:
    return getPtEnumeratorSingle(helm, params.placeholders);
  case PolyToolEnumeratorTypes.Matrix:
    return getPtEnumeratorMatrix(helm, params.placeholders);
  }
}
