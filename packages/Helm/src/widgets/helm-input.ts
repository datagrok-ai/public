import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {fromEvent, Observable, Unsubscribable} from 'rxjs';

import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmInputBase, IHelmHelper, IHelmInputInitOptions} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmAtom, HelmMol, IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {SeqValueBase} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {MonomerNumberingTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {defaultErrorHandler} from '../utils/err-info';
import {getHoveredMonomerFromEditorMol, getSeqMonomerFromHelmAtom} from '../utils/get-hovered';


import {_package} from '../package';

/**
 * Default editor size (px). The editor canvas AND the host element both default
 * to this so `ui.input.helmAsync` renders a compact, predictable box instead of
 * stretching to the (often wide) form row or falling back to hwe's 400×100
 * editor default. Override per-input via `editorOptions.{width,height}`, or — for
 * a fill layout — set `width/height` directly on `getInput()` (e.g. `'100%'`),
 * which wins over these defaults; `onSizeChanged` then re-fits the canvas.
 */
const DEFAULT_HELM_INPUT_HEIGHT = 250;
const DEFAULT_HELM_INPUT_WIDTH = 250;

export class HelmInput extends HelmInputBase {
  /** Input type identifier (such as "Slider" for the slider input). See {@link InputType}. */
  get inputType(): string {
    return 'Macromolecule';
  }

  /** Data type this input can edit. See {@link DG.Type}. */
  get dataType(): string {
    return DG.TYPE.STRING;
  }

  getInput(): HTMLElement {
    return this.viewerHost;
  }

  private _seqValue: SeqValueBase;

  getValue(): SeqValueBase { return this._seqValue; }

  setValue(value: SeqValueBase): void {
    this._seqValue = value;
    this.viewer.editor.setHelm(this._seqValue.helm);
  }

  getStringValue(): string { return this._seqValue?.helm ?? ''; }

  setStringValue(value: string): void {
    if (!this._seqValue) {
      let seqCol: DG.Column;
      DG.DataFrame.fromColumns([seqCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [value])]);
      seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
      seqCol.meta.units = NOTATION.HELM;
      const sh = this.helmHelper.seqHelper.getSeqHandler(seqCol);
      this._seqValue = sh.getValue(0);
    }
    this.viewer.editor.setHelm(this._seqValue.value = value);
  }

  // -- IHelmInput --

  get molValue(): HelmMol {
    return this.viewer.editor.m;
  }

  set molValue(value: HelmMol) {
    this.viewer.editor.setMol(value);
    this._seqValue.value = this.viewer.editor.getHelm();
  }

  get onMouseMove(): Observable<MouseEvent> { return fromEvent<MouseEvent>(this.viewer.host, 'mousemove'); }

  get onClick(): Observable<MouseEvent> { return fromEvent<MouseEvent>(this.viewer.host, 'click'); }

  readonly viewerHost: HTMLDivElement;
  readonly viewer: IHelmWebEditor;
  readonly editHintDiv: HTMLDivElement;

  helmString: string = '';
  helmSelection: number[];

  private subs: Unsubscribable[];

  protected constructor(
    private readonly helmHelper: IHelmHelper,
    private readonly monomerLib: IMonomerLibBase,
    name?: string,
    public readonly options?: IHelmInputInitOptions,
    private readonly logger = _package.logger,
  ) {
    super();

    if (name) this.captionLabel.innerText = name;

    // The host defaults to a fixed DEFAULT_HELM_INPUT_{WIDTH,HEIGHT} box (NOT
    // `width:100%`, which would stretch to the wide form row and surface the
    // hwe 400×100 editor default). A consumer that wants a fill layout sets
    // `width/height` on `getInput()` directly (those inline overrides win); the
    // `min-*` floor keeps it visible if they set `width:100%` in a 0-size
    // parent. The editor canvas is created at the SAME default so the very first
    // paint is already the right size (no 400-wide flash), and `onSizeChanged`
    // below re-fits it from the host's own box whenever the host is resized.
    const defaultHeight = options?.editorOptions?.height ?? DEFAULT_HELM_INPUT_HEIGHT;
    const defaultWidth = options?.editorOptions?.width ?? DEFAULT_HELM_INPUT_WIDTH;
    this.viewerHost = ui.div([], {
      classes: 'ui-input-editor',
      style: {
        width: `${defaultWidth}px`,
        height: `${defaultHeight}px`,
        minWidth: `${defaultWidth}px`,
        minHeight: `${defaultHeight}px`,
        overflow: 'hidden',
      },
    });
    this.viewer = this.helmHelper.createHelmWebEditor(this.viewerHost, {
      ...options?.editorOptions,
      width: defaultWidth,
      height: defaultHeight,
    });

    this.editHintDiv = ui.divH([
      ui.link('Click to edit', () => { this.viewer.host.click(); }, undefined, {}),
    ], {style: {display: 'none', position: 'absolute', top: '0', right: '0'}});

    this.subs = [];
    this.subs.push(this.monomerLib.onChanged
      .subscribe(() => this.viewer.editor.redraw()));
    /* eslint-disable rxjs/no-ignored-subscription */
    fromEvent<MouseEvent>(this.viewer.host, 'mousemove')
      .subscribe(this.viewerOnMouseMove.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'mouseenter')
      .subscribe(this.viewerOnMouseEnter.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'mouseleave')
      .subscribe(this.viewerOnMouseLeave.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'click')
      .subscribe(this.viewerOnClick.bind(this));
    /* eslint-enable rxjs/no-ignored-subscription */

    // Drive the editor canvas from the host element's OWN measured box so it
    // works no matter how the consumer mounts the input — via `root` (form
    // layout) or by lifting `getInput()` (the host) into a custom container.
    // (The old code read `this.root.clientHeight`, which is 0 whenever only
    // `getInput()` is in the DOM — producing the reported 0-height editor.)
    this.subs.push(ui.onSizeChanged(this.input).subscribe(() => {
      const w = this.input.clientWidth;
      const h = this.input.clientHeight;
      if (w > 0 && h > 0)
        this.viewer.editor.setSize(w, h);
    }));

    this.root.classList.add('ui-input-helm');
    this.root.append(this.viewerHost, this.editHintDiv);
  }

  detach(): void {
    for (const sub of this.subs)
      sub.unsubscribe();
  }

  protected toLog(): string {
    return `Helm: HelmInput<#>`;
  }

  /** Inspired by {@link DG.MarkdownInput}, {@link ui._create} */
  static create(helmHelper: IHelmHelper, monomerLib: IMonomerLibBase,
    name?: string, options?: IHelmInputInitOptions
  ): HelmInput {
    const input = new HelmInput(helmHelper, monomerLib, name, options);
    // TODO: Apply options
    if (options?.value != null)
      input.value = options.value;
    return input;
  }

  redraw(): void {
    this.viewer.editor.redraw();
  }

  showTooltip(content: HTMLElement | string, a: HelmAtom): void {
    /** Monomer sign offset for tooltip inspired by {@link org.helm.webeditor.drawMonomer()} */
    const so = this.viewer.editor.getDrawOptions().fontsize * 1.2;
    const bcr = this.input.getBoundingClientRect();
    ui.tooltip.show(content, bcr.left + a.p.x + so, bcr.top + a.p.y + so);
  }

  // -- Handle events --

  private viewerOnClick(_event: MouseEvent): void {
    ui.tooltip.hide();

    if (this.options?.editable != false)
      this.showEditorDialog();
  }

  private viewerOnMouseMove(event: MouseEvent): void {
    const logPrefix = `${this.toLog()}.viewerOnMouseMove()`;

    const argsX = event.offsetX;
    const argsY = event.offsetY;
    const mol = this.viewer.editor.m;
    const hoveredAtom = getHoveredMonomerFromEditorMol(argsX, argsY, mol, this.viewer.editor.div.clientHeight);
    if (hoveredAtom) {
      const seqMonomer = getSeqMonomerFromHelmAtom(this.value, hoveredAtom);
      const monomerLib = _package.monomerLib;
      const tooltipEl = monomerLib ? monomerLib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) :
        ui.divText('Monomer library is not available');
      this.showTooltip(tooltipEl, hoveredAtom);
      event.preventDefault();
      event.stopPropagation();
    } else
      ui.tooltip.hide();

    //this.logger.debug(`${logPrefix}, x = ${event.x}, y = ${event.y}`);
  }

  private viewerOnMouseEnter(_event: MouseEvent): void {
    //this.logger.debug(`${this.toLog()}.viewerOnMouseEnter()`);
    if (this.options?.editable != false)
      this.showEditHint();
  }

  private viewerOnMouseLeave(event: MouseEvent): void {
    if (this.editHintDiv.contains(event.relatedTarget as Node)) return; // prevents flickering

    //this.logger.debug(`${this.toLog()}.viewerOnMouseLeave()`);
    if (this.options?.editable != false)
      this.hideEditHint();
  }

  // editHint

  showEditHint(): void {
    this.editHintDiv.style.removeProperty('display');
  }

  hideEditHint(): void {
    this.editHintDiv.style.display = 'none';
  }

  // editor dialog

  showEditorDialog(): void {
    const webEditorHost = ui.div();
    webEditorHost.style.height = '100%';
    const webEditorApp = this.helmHelper.createWebEditorApp(webEditorHost, this.viewer.editor.getHelm());

    const destroyEditorDialog = () => {
      $(webEditorHost).empty();
    };

    // Edit the molecule and select monomers to enumerate
    const _editorDialog = ui.dialog({showHeader: false, showFooter: true})
      .add(webEditorHost!)
      .onOK(() => {
        try {
          const webEditorValue = webEditorApp!.canvas!.getHelm(true)
            .replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
          this.viewer.editor.setHelm(webEditorValue);
          this.helmString = webEditorValue;

          const selection = webEditorApp!.canvas!.helm!.jsd.m.atoms;
          this.helmSelection = [];
          for (let i = 0; i < selection.length; i++) {
            if (selection[i].selected)
              this.helmSelection.push(i);

            this.viewer.editor.m.atoms[i].highlighted = selection[i].selected;
          }
          this.viewer.editor.redraw();
        } catch (err: any) {
          defaultErrorHandler(err);
        } finally {
          destroyEditorDialog();
        }
      })
      .onCancel(() => {
        destroyEditorDialog();
      })
      .show({modal: true, fullScreen: true});
  }
}
