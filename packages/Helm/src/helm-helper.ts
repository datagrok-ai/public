import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as scil from 'scil';
import * as org from 'org';

import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';

import {HelmWebEditor} from './helm-web-editor';

type HelmHelperWindowType = Window & {
  $helmHelper?: HelmHelper,
}
declare const window: HelmHelperWindowType;

export class HelmHelper implements IHelmHelper {
  createHelmWebEditor(): IHelmWebEditor {
    return new HelmWebEditor();
  }

  createWebEditorApp(host: HTMLDivElement, helm: string): org.helm.IWebEditorApp {
    org.helm.webeditor.MolViewer.molscale = 0.8;
    const webEditorApp: org.helm.IWebEditorApp = new org.helm.webeditor.App(host, {
      showabout: false,
      mexfontsize: '90%',
      mexrnapinontab: true,
      topmargin: 20,
      mexmonomerstab: true,
      sequenceviewonly: false,
      mexfavoritefirst: true,
      mexfilter: true
    });
    const sizes = webEditorApp.calculateSizes();
    webEditorApp.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
    let s: Partial<CSSStyleDeclaration> = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
    scil.apply(webEditorApp.sequence.style, s);
    scil.apply(webEditorApp.notation.style, s);
    s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + webEditorApp.toolbarheight) + 'px'};
    scil.apply(webEditorApp.properties.parent.style, s);
    webEditorApp.structureview.resize(sizes.rightwidth, sizes.bottomheight + webEditorApp.toolbarheight);
    webEditorApp.mex.resize(sizes.topheight - 80);
    setTimeout(() => {
      webEditorApp.canvas.helm.setSequence(helm, 'HELM');
    }, 200);
    return webEditorApp;
  }

  // -- Instance singleton --

  public static async getInstance(): Promise<IHelmHelper> {
    let res: HelmHelper = window.$helmHelper!;
    if (!res) {
      // create & init singleton
      window.$helmHelper = res = new HelmHelper();
      // await res.init();
    }
    return res;
  }
}
