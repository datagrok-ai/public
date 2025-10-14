import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {Atom, HelmType, IWebEditorMonomer, GetMonomerResType, IHelmBio, HelmAtom} from '@datagrok-libraries/bio/src/helm/types';

import {OrgHelmModule, ScilModule} from '../types';
import {RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from '../constants';

import {_package} from '../package';

declare const org: OrgHelmModule;
declare const scil: ScilModule;

type GetMonomerOverridingFunc = (
  a: HelmAtom | HelmType, name: string | undefined, monomerLib: IMonomerLib,
  originalGetMonomer: (a: HelmAtom | HelmType, name: string) => GetMonomerResType) => GetMonomerResType;

export function getMonomerOverrideAndLogAlert(
  monomerLib: IMonomerLib, getMonomerOverriding: GetMonomerOverridingFunc, trigger: () => void, logger: ILogger,
): void {
  const logPrefix = 'Helm: overrideGetMonomerAndAlert()';
  const monomers = org.helm.webeditor.Monomers;
  const getMonomerOriginal = monomers.getMonomer.bind(monomers);
  const alertOriginal = scil.Utils.alert;
  try {
    org.helm.webeditor.Monomers.getMonomer =
      (a: HelmAtom | HelmType, name?: string): GetMonomerResType => {
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
        const rs = monomerLib.getRS(monomer[SMILES].toString());
        if(rs == null || Object.keys(rs).length === 0) {
          isBroken = true;
        } else {
          webEditorMonomer.rs = Object.keys(rs).length;
          webEditorMonomer.at = rs;
        }
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
