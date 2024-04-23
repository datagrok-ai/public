import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {Subject, Unsubscribable} from 'rxjs';

import {HelmProps, HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {svgToImage} from '@datagrok-libraries/utils/src/svg';

import {IEditor} from '../helm-monomer-placer';

import {_package} from '../package';

export class HelmService extends HelmServiceBase {
  private readonly hostDiv: HTMLDivElement;

  private editor: IEditor | null = null;
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

  protected override async requestRender(
    key: keyof any | undefined, task: RenderTask<HelmProps>, renderHandler: () => void
  ): Promise<[Unsubscribable, number, () => void]> {
    const logPrefix = `${this.toLog()}.requestRender()`;
    this.logger.debug(`${logPrefix}, ` + `key: ${key?.toString()}`);
    const emptyCanvasHash: number = 0;

    if (!this.editor) {
      this.editor = new JSDraw2.Editor(this.hostDiv, {width: 300, height: 300, skin: 'w8', viewonly: true}) as IEditor;
    }
    const lST = window.performance.now();
    this.editor.resize(task.props.width, task.props.height);
    this.editor.setData(task.props.gridCell.cell.valueString, 'helm');
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
      renderHandler();
    };

    const renderEvent = new Subject();
    const renderSub = renderEvent.subscribe(() => {
      const logPrefixInt = `${logPrefix} on renderEvent`;
      this.logger.debug(`${logPrefixInt}`);
    });
    const trigger = (): void => {
      svgToImage(svgEl, dpr).then((imageEl) => {
        this.image = imageEl;
        renderHandlerInt();
      });
    };
    return [renderSub, -1, trigger];
  }

  protected onRendered(
    key: keyof any | undefined, task: RenderTask<HelmProps>, emptyCanvasHash: number
  ): boolean {
    const dpr = window.devicePixelRatio;
    const canvas = ui.canvas(task.props.width, task.props.height);
    try {
      const g = canvas.getContext('2d');
      if (!this.image || !g) return false;
      // g.fillStyle = '#C0C0FF';
      // g.fillRect(0, 0, canvas.width, canvas.height);
      g.drawImage(this.image, 0, 0);
      task.onAfterRender(canvas);
    } finally {
      canvas.remove();
    }
    return true;
  }

  async reset(): Promise<void> { }
}
