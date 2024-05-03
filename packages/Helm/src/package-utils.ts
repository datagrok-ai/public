import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import org from 'org';
import scil from 'scil';

import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {HelmService} from './utils/helm-service';
import {GetMonomerFunc, GetMonomerResType} from './utils/get-monomer';

type HelmWindowType = Window & {
  $helmService?: HelmServiceBase,
}
declare const window: HelmWindowType;

export function _getHelmService(): HelmServiceBase {
  let res = window.$helmService;
  if (!res) res = window.$helmService = new HelmService();
  return res;
}

export function initHelmPatchDojo(): void {
  // patch window.dojox.gfx.svg.Text.prototype.getTextWidth hangs
  /** get the text width in pixels */
  // @ts-ignore
  window.dojox.gfx.svg.Text.prototype.getTextWidth = function() {
    const rawNode = this.rawNode;
    const oldParent = rawNode.parentNode;
    const _measurementNode = rawNode.cloneNode(true);
    _measurementNode.style.visibility = 'hidden';

    // solution to the "orphan issue" in FF
    let _width = 0;
    const _text = _measurementNode.firstChild.nodeValue;
    oldParent.appendChild(_measurementNode);

    // solution to the "orphan issue" in Opera
    // (nodeValue == "" hangs firefox)
    if (_text != '') {
      let watchdogCounter = 100;
      while (!_width && --watchdogCounter > 0) { // <-- hangs
        //Yang: work around svgweb bug 417 -- http://code.google.com/p/svgweb/issues/detail?id=417
        if (_measurementNode.getBBox)
          _width = parseInt(_measurementNode.getBBox().width);
        else
          _width = 68;
      }
    }
    oldParent.removeChild(_measurementNode);
    return _width;
  };
}

/** Patches Pistoia Monomers.getMonomer method to utilize DG Bio monomer Lib */
export function initHelmPatchPistoia(monomerLib: IMonomerLib, logger: ILogger): void {
  const logPrefix: string = 'Helm: initHelmPatchPistoia()';
  const monomers = org.helm.webeditor.Monomers;
  const getMonomerOriginal: GetMonomerFunc = monomers.getMonomer.bind(monomers);
  const alertOriginal = scil.Utils.alert;
  org.helm.webeditor.Monomers.getMonomer = (
    a: org.helm.IAtom | org.helm.HelmType, name: string
  ): GetMonomerResType => {
    const oWem = getMonomerOriginal(a, name);
    const dgWem = monomerLib.getWebEditorMonomer(a, name);
    return oWem;
  };
  scil.Utils.alert = (s: string): void => {
    logger.warning(`${logPrefix}, scil.Utils.alert() s = 's'.`);
  };

  org.helm.webeditor.Monomers = new class {
    constructor(base: org.helm.IMonomers) {
      return new Proxy(base, {
        get(target: any, p: string | symbol, _receiver: any): any {
          return target[p];
        },
        set(target: any, p: string | symbol, newValue: any, _receiver: any): boolean {
          target[p] = newValue;
          return true;
        },
        apply(target: any, thisArg: any, argArray: any[]): any {
          return target[thisArg](...argArray);
        },
      });
    }
  }(org.helm.webeditor.Monomers) as org.helm.IMonomers;

  org.helm.webeditor = new class {
    constructor(base: org.helm.IOrgHelmWebEditor) {
      return new Proxy(base, {
        get(target: any, p: string | symbol, _receiver: any): any {
          return target[p];
        },
        set(target: any, p: string | symbol, newValue: any, _receiver: any): boolean {
          target[p] = newValue;
          return true;
        },
        apply(target: any, thisArg: any, argArray: any[]): any {
          return target[thisArg](...argArray);
        },
      });
    }
  }(org.helm.webeditor) as org.helm.IOrgHelmWebEditor;
}
