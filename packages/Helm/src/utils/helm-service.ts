import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {Subject, Unsubscribable} from 'rxjs';
import {LRUCache} from 'lru-cache';

import {HelmEditor, HelmType, IHelmBio, IHelmEditorOptions} from '@datagrok-libraries/bio/src/helm/types';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {HelmAux, HelmProps, HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {svgToImage} from '@datagrok-libraries/utils/src/svg';

import {JSDraw2Module} from '../types';

import {_package} from '../package';

declare const JSDraw2: JSDraw2Module;

export class HelmService extends HelmServiceBase {
  private readonly hostDiv: HTMLDivElement;

  private readonly editorLruCache: LRUCache<string, HelmEditor>;

  // private editor!: HelmEditor;
  private image: HTMLImageElement | null = null;

  constructor() {
    super(_package.logger);

    this.hostDiv = ui.box();
    // this.hostDiv.style.display = 'none'; // Disables drawing at all
    this.hostDiv.style.position = 'fixed';
    this.hostDiv.style.left = '0px';
    this.hostDiv.style.right = '0px';
    this.hostDiv.style.width = '0px';
    this.hostDiv.style.height = '0px';
    this.hostDiv.style.visibility = 'hidden';
    document.body.appendChild(this.hostDiv);

    this.editorLruCache = new LRUCache<string, HelmEditor>({max: 20});
  }

  protected override toLog(): string {
    return `Helm: ${super.toLog()}`;
  }

  private getEditor(props: HelmProps): HelmEditor {
    const editorKey = props.monomerLib.source;

    let resEditor: HelmEditor | undefined = this.editorLruCache.get(editorKey);
    if (!resEditor) {
      const getMonomerFuncs = _package.helmHelper.buildMonomersFuncsFromLib(props.monomerLib, );
      resEditor = new JSDraw2.Editor<HelmType, IHelmBio, IHelmEditorOptions>(this.hostDiv, {
        width: props.width, height: props.height, skin: 'w8', viewonly: true,
        drawOptions: {getMonomer: getMonomerFuncs.getMonomer},
      });
      this.editorLruCache.set(editorKey, resEditor);
    }
    return resEditor;
  }

  protected override async requestRender(
    key: keyof any | undefined, task: RenderTask<HelmProps, HelmAux>, renderHandler: (aux: HelmAux) => void
  ): Promise<[Unsubscribable, number, () => void]> {
    const logPrefix = `${this.toLog()}.requestRender()`;
    this.logger.debug(`${logPrefix}, start, ` + `key: ${key?.toString()}`);
    const emptyCanvasHash: number = 0;
    const monomerLib = task.props.monomerLib;

    const editor = this.getEditor(task.props);

    const lST = window.performance.now();
    editor.options.width = task.props.width;
    editor.options.height = task.props.height;
    editor.resize(task.props.width, task.props.height);
    const helmStr = task.props.helm;
    if (helmStr) {
      // getMonomerOverrideAndLogAlert(
      //   this.monomerLib,
      //   getMonomerPatched,
      //   () => {
      //     this.logger.debug(`${logPrefix}, editor.setData( '${helmStr}' )`);
      //     this.editor.setData(helmStr, 'helm');
      //   }, this.logger);
      this.logger.debug(`${logPrefix}, editor.setData( '${helmStr}' )`);
      editor.setData(helmStr, 'helm');
    }
    if (!helmStr)
      editor.reset();

    const svgEl = editor.div.children[0] as SVGSVGElement;
    const bBox = svgEl.getBBox();
    const [cellWidth, cellHeight] = [task.props.width, task.props.height];
    const bScale = Math.min(cellWidth * 0.95 / bBox.width, cellHeight * 0.95 / bBox.height);
    const [drawWidth, drawHeight] = [bBox.width * bScale, bBox.height * bScale];
    const aux: HelmAux = {
      mol: editor.m.clone(false),
      bBox: new DG.Rect(bBox.x, bBox.y, bBox.width + 2, bBox.height + 2) /* adjust for clipping right/bottom */,
      dBox: new DG.Rect((cellWidth - drawWidth) / 2, (cellHeight - drawHeight) / 2, drawWidth, drawHeight),
      cBox: new DG.Rect(0, 0, task.props.width, task.props.height),
    };
    const lET = window.performance.now();

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
      g.drawImage(this.image,
        aux.bBox.x, aux.bBox.y, aux.bBox.width, aux.bBox.height,
        aux.dBox.x, aux.dBox.y, aux.dBox.width, aux.dBox.height);
      task.onAfterRender(canvas, aux);
    } finally {
      canvas.remove();
    }
    return true;
  }

  async reset(): Promise<void> { }
}

