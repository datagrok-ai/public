import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as org from 'org';
import * as scil from 'scil';
import * as JSDraw2 from 'JSDraw2';

import $ from 'cash-dom';
import {Subject, Unsubscribable} from 'rxjs';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {HelmAux, HelmProps, HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {svgToImage} from '@datagrok-libraries/utils/src/svg';

import {_package} from '../package';
import {DummyWebEditorMonomer} from './dummy-monomer';

export class HelmService extends HelmServiceBase {
  private readonly hostDiv: HTMLDivElement;

  private editor: JSDraw2.Editor | null = null;
  private image: HTMLImageElement | null = null;

  constructor() {
    super(_package.logger);

    this.hostDiv = ui.box();
    // this.hostDiv.style.display = 'none'; // Disables drawing at all
    this.hostDiv.style.position = 'absolute';
    this.hostDiv.style.left = '0px';
    this.hostDiv.style.right = '0px';
    this.hostDiv.style.width = '0px';
    this.hostDiv.style.height = '0px';
    this.hostDiv.style.visibility = 'hidden';
    document.body.appendChild(this.hostDiv);
  }

  protected override toLog(): string {
    return `Helm: ${super.toLog()}`;
  }

  protected override async requestRender(
    key: keyof any | undefined, task: RenderTask<HelmProps, HelmAux>, renderHandler: (aux: HelmAux) => void
  ): Promise<[Unsubscribable, number, () => void]> {
    const logPrefix = `${this.toLog()}.requestRender()`;
    this.logger.debug(`${logPrefix}, start, ` + `key: ${key?.toString()}`);
    const emptyCanvasHash: number = 0;

    if (!this.editor) {
      this.editor = new JSDraw2.Editor(this.hostDiv,
        {width: task.props.width, height: task.props.height, skin: 'w8', viewonly: true});
    }
    const lST = window.performance.now();
    this.editor.options.width = task.props.width;
    this.editor.options.height = task.props.height;
    this.editor.resize(task.props.width, task.props.height);
    const helmStr = task.props.helm;
    let hasMissing: boolean = false;
    if (helmStr) {
      const monomers = org.helm.webeditor.Monomers;
      const originalGetMonomer = monomers.getMonomer;
      const originalAlert = scil.Utils.alert;
      try {
        org.helm.webeditor.Monomers.getMonomer = (a: org.helm.IAtom | string, name: string): org.helm.IAtom | null => {
          let resAtom: any = originalGetMonomer.bind(monomers)(a, name);

          if (!resAtom) {
            // Input logic
            if (a == null && name == null)
              return null;

            let s: string;
            let biotype: string;
            if (name == null) {
              biotype = (a as org.helm.IAtom).biotype();
              s = (a as org.helm.IAtom).elem;
            } else {
              biotype = a as string;
              s = org.helm.webeditor.IO.trimBracket(name);
            }
            resAtom = new DummyWebEditorMonomer(biotype, s);
            //hasMissing = true;
          }
          return resAtom;
        };
        // Preventing alert message box for missing monomers with compressed Scilligence.JSDraw2.Lite.js
        scil.Utils.alert = (s: string): void => {
          this.logger.warning(`${logPrefix}, scil.Utils.alert() s = 's'.`);
        };
        this.logger.debug(`${logPrefix}, editor.setData( '${helmStr}' )`);
        this.editor.setData(helmStr, 'helm');
      } finally {
        monomers.getMonomer = originalGetMonomer;
        scil.Utils.alert = originalAlert;
      }
    }
    if (!helmStr || hasMissing)
      this.editor.reset();

    const bBox = (this.editor.div.children[0] as SVGSVGElement).getBBox();
    const aux: HelmAux = {
      mol: this.editor.m.clone(false),
      bBox: new DG.Rect(bBox.x, bBox.y, bBox.width + 2, bBox.height + 2) /* adjust for clipping right/bottom */,
      cBox: new DG.Rect(0, 0, task.props.width, task.props.height),
    };
    const lET = window.performance.now();

    const svgEl = $(this.hostDiv).find('svg').get(0) as unknown as SVGSVGElement;
    const dpr = window.devicePixelRatio;

    const rST = window.performance.now();
    const renderHandlerInt = (): void => {
      const rET: number = window.performance.now();
      this.logger.debug(`${logPrefix}.renderHandlerInt(), ` +
        `key: ${key?.toString()}, ` +
        `load: ${lET - lST} ms,\n    ` + `renderET: ${rET - rST} ms, ` +
        `emptyCanvasHash: ${emptyCanvasHash}`);
      renderHandler(aux);
    };

    const renderEvent = new Subject();
    const renderSub = renderEvent.subscribe(() => {
      const logPrefixInt = `${logPrefix} on renderEvent`;
      this.logger.debug(`${logPrefixInt}`);
    });
    const trigger = (): void => {
      svgToImage(svgEl, dpr).then((imageEl) => {
        this.image = imageEl;
        renderHandlerInt(); // calls onRendered()
      });
    };
    this.logger.debug(`${logPrefix}, end, ` + `key: ${key?.toString()}`);
    return [renderSub, -1, trigger];
  }

  protected onRendered(
    key: keyof any | undefined, task: RenderTask<HelmProps, HelmAux>, emptyCanvasHash: number, aux: HelmAux
  ): boolean {
    const dpr = window.devicePixelRatio;
    const canvas = ui.canvas(task.props.width, task.props.height);
    try {
      const g = canvas.getContext('2d');
      if (!this.image || !g) return false;
      // g.fillStyle = '#C0C0FF';
      // g.fillRect(0, 0, canvas.width, canvas.height);
      g.drawImage(this.image, 0, 0);
      task.onAfterRender(canvas, aux);
    } finally {
      canvas.remove();
    }
    return true;
  }

  async reset(): Promise<void> { }
}
