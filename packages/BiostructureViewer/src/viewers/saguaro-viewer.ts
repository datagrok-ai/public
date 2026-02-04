import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {Observable, Subject, Unsubscribable} from 'rxjs';
import {v4 as uuidv4} from 'uuid';
import {RcsbFv, RcsbFvDisplayTypes, RcsbFvRowConfigInterface} from '@rcsb/rcsb-saguaro';
import {RcsbFvBoardConfigInterface} from '@rcsb/rcsb-saguaro/build/RcsbFv/RcsbFvConfig/RcsbFvConfigInterface';

import {intToHtmlA} from '@datagrok-libraries/utils/src/color';
import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {testEvent} from '@datagrok-libraries/test/src/test';
import {BiotrackProps} from '@datagrok-libraries/bio/src/viewers/biotrack';

import {_package} from '../package';

const enum PROPS_CATS {
  DATA = 'Data',
  STYLE = 'Style',
  LAYOUT = 'Layout',
  CONTROLS = 'Controls',
}

export enum PROPS {
  // -- Data --
  // -- Style --
  // -- Layout --
  rowTitleWidth = 'rowTitleWidth',
  trackWidth = 'trackWidth',
  includeAxis = 'includeAxis',
  includeTooltip = 'includeTooltip',
  disableMenu = 'disableMenu',
  borderColor = 'borderColor',
  borderWidth = 'borderWidth',
  hideInnerBorder = 'hideInnerBorder',
  hideTrackFrameGlow = 'hideTrackFrameGlow',
  highlightHoverPosition = 'highlightHoverPosition',
  hideRowGlow = 'hideRowGlow',

  // -- Controls --
}

/** rcsb molstar track viewer */
export class SaguaroViewer extends DG.JsViewer {
  private viewed: boolean = false;

  // -- Data --
  // -- Style --
  // -- Layout --
  [PROPS.rowTitleWidth]?: number;
  [PROPS.trackWidth]: number;
  [PROPS.includeAxis]: boolean;
  [PROPS.includeTooltip]: boolean;
  [PROPS.disableMenu]: boolean;
  [PROPS.borderColor]: number;
  [PROPS.borderWidth]: number;
  [PROPS.hideInnerBorder]: boolean;
  [PROPS.hideTrackFrameGlow]: boolean;
  [PROPS.highlightHoverPosition]: boolean;
  [PROPS.hideRowGlow]: boolean;

  // -- Controls --

  constructor() {
    super();

    // -- Data --
    // -- Style --
    // -- Layout --
    this.rowTitleWidth = this.int(PROPS.rowTitleWidth, 75,
      {category: PROPS_CATS.LAYOUT, editor: 'slider', min: 0, max: 75});
    this.trackWidth = this.int(PROPS.trackWidth, 16,
      {category: PROPS_CATS.LAYOUT, editor: 'slider', min: 4, max: 200});

    // -- Style --
    this.borderColor = this.int(PROPS.borderColor, DG.Color.lightGray,
      {category: PROPS_CATS.STYLE});
    this.borderWidth = this.int(PROPS.borderWidth, 0,
      {category: PROPS_CATS.STYLE, editor: 'slider', min: 0, max: 16});
    this.hideInnerBorder = this.bool(PROPS.hideInnerBorder, true,
      {category: PROPS_CATS.STYLE});
    this.hideTrackFrameGlow = this.bool(PROPS.hideTrackFrameGlow, false,
      {category: PROPS_CATS.STYLE});
    this.highlightHoverPosition = this.bool(PROPS.highlightHoverPosition, true,
      {category: PROPS_CATS.STYLE});

    // -- Controls --
    this.includeTooltip = this.bool(PROPS.includeTooltip, false,
      {category: PROPS_CATS.LAYOUT});
    this.disableMenu = this.bool(PROPS.disableMenu, false,
      {category: PROPS_CATS.LAYOUT});

    this.subs.push(
      ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    this.viewSyncer = new PromiseSyncer(_package.logger);
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = ++SaguaroViewer.viewerCounter;

  private viewerToLog(): string { return `MolstarViewer<${this.viewerId}>`; }

  override onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (!property) {
      _package.logger.warning('SaguaroViewer.onPropertyChanged() property is null');
      return;
    }
    switch (property.name) {
      case PROPS.borderColor:
        const _borderColorValue: string = intToHtmlA(this.borderColor);
        break;

      case PROPS.rowTitleWidth:
      case PROPS.trackWidth:
      case PROPS.includeAxis:
      case PROPS.includeTooltip:
      case PROPS.disableMenu:
      case PROPS.borderWidth:
      case PROPS.hideInnerBorder:
      case PROPS.hideTrackFrameGlow:
      case PROPS.highlightHoverPosition:
      case PROPS.hideRowGlow:
        const propValue: any = this.props.get(property.name);
        this.boardConfigDataProps[property.name as keyof RcsbFvBoardConfigInterface] = propValue;
        break;
    }
  }

  override onTableAttached(): void {
    _package.logger.debug('SaguaroViewer.onTableAttached(), ');
    const superOnTableAttached = super.onTableAttached.bind(this);

    // -- Props editors --

    superOnTableAttached();
    this.setData();
  }

  override detach(): void {
    const superDetach = super.detach.bind(this);
    this.viewSyncer.sync('detach()', async () => {
      if (this.setDataInProgress) return; // check setDataInProgress synced
      if (this.viewed) {
        await this.destroyView('detach'); // detach
        this.viewed = false;
      }
      superDetach();
    });
  }

  // -- Data --

  setData(): void {
    this.viewSyncer.sync('setData()', async () => { // setData
      if (!this.setDataInProgress) this.setDataInProgress = true; else return; // check setDataInProgress synced
      try {
        if (this.viewed) {
          await this.destroyView('setData'); // setData
          this.viewed = false;
        }
        // TODO: Data

        if (!this.viewed) {
          await ui.tools.waitForElementInDom(this.root);
          await this.buildView('setData'); // setData
          this.viewed = true;
        }
      } finally {
        this.setDataInProgress = false;
      }
    });
  }

  // -- View --
  private viewSyncer;
  private setDataInProgress: boolean = false;

  private viewerDivId: string;
  private viewerDiv?: HTMLDivElement;
  private viewer?: RcsbFv;
  private splashDiv?: HTMLDivElement;
  private viewSubs: Unsubscribable[] = [];

  // properties containers
  boardConfigDataProps: RcsbFvBoardConfigInterface = {};

  private async destroyView(purpose: string): Promise<void> {
    _package.logger.debug(`BiotrackViewer.destroyView( purpose='${purpose}' )`);
    if (true) {
      // Clear viewer
    }

    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }
  }

  private async buildView(purpose: string): Promise<void> {
    _package.logger.debug(`BiotrackViewer.buildView( purpose='${purpose}' )`);

    if (!this.viewerDiv) {
      this.viewerDivId = uuidv4();
      this.viewerDiv = ui.div([], {
        id: this.viewerDivId,
        classes: 'd4-saguaro-viewer',
        style: {width: '100%', height: '100%', backgroundColor: '#FAFAFF'},
      });
      this.root.appendChild(this.viewerDiv);

      const sequence =
        'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQ' +
        'EEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA' +
        'RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMS';


      const boardConfig = {
        length: 110,
        range: {
          min: 1,
          max: sequence.length,
        },
        includeAxis: true,
      };

      const sequenceTrack: RcsbFvRowConfigInterface = {
        trackId: 'sequenceTrack',
        trackHeight: 20,
        trackColor: '#F9F9F9',
        displayType: RcsbFvDisplayTypes.SEQUENCE,
        nonEmptyDisplay: true,
        rowTitle: 'SEQUENCE',
        trackData: [
          {
            begin: 1,
            value: sequence,
          },
        ],
      };

      const bondTrack: RcsbFvRowConfigInterface = {
        trackId: 'bondTrack',
        trackHeight: 20,
        trackColor: '#F9F9F9',
        displayType: RcsbFvDisplayTypes.BOND,
        displayColor: '#6572E3',
        rowTitle: 'BOND',
        trackData: [
          {
            begin: 30,
            end: 90,
            isEmpty: true,
          },
        ],
      };

      const blockTrack: RcsbFvRowConfigInterface = {
        trackId: 'blockTrack',
        trackHeight: 20,
        trackColor: '#F9F9F9',
        displayType: RcsbFvDisplayTypes.BLOCK,
        displayColor: '#FF0000',
        rowTitle: 'BLOCK',
        trackData: [
          {
            begin: 30,
            end: 60,
            gaps: [{
              begin: 40,
              end: 50,
              isConnected: true,
            }],
          }, {
            begin: 80,
            end: 90,
            openEnd: true,
          },
        ],
      };

      const pinTrack: RcsbFvRowConfigInterface = {

        trackId: 'pinTrack',
        trackHeight: 20,
        trackColor: '#F9F9F9',
        displayType: RcsbFvDisplayTypes.PIN,
        displayColor: '#65F253',
        rowTitle: 'PIN',
        trackData: [
          {
            begin: 60,
          }, {
            begin: 70,
          }, {
            begin: 100,
          },
        ],
      };

      const vline: RcsbFvRowConfigInterface = {
        trackId: 'vlineTrack',
        trackHeight: 20,
        trackColor: '#F9F9F9',
        displayType: RcsbFvDisplayTypes.VLINE,
        displayColor: '#65F253',
        rowTitle: 'VLINE',
      };

      const rowConfig: Array<RcsbFvRowConfigInterface> = [sequenceTrack, bondTrack, blockTrack, pinTrack, vline];

      this.viewer = new RcsbFv({
        boardConfigData: boardConfig,
        rowConfigData: rowConfig,
        elementId: this.viewerDivId,
      });
    }
  }

  private calcSize(): void {
    this.viewSyncer.sync('calcSize()', async () => {
      if (!this.viewer || !this.viewerDiv) return;

      const cw: number = this.root.clientWidth;
      const ch: number = this.root.clientHeight;
      _package.logger.debug('BiotrackViewer.calcSize( ${cw.toString()} x ${ch.toString()} )');

      this.viewerDiv.style.width = `${cw}px`;
      this.viewerDiv.style.height = `${ch}px`;

      // viewer handle resize
      await this.viewer.updateBoardConfig({
        boardConfigData:
          {
            trackWidth: cw - (this.rowTitleWidth ?? 0),
            includeAxis: true,
          },
      }).then(() => {
        // -- Props editors --
        const rowTitleWidthProp: DG.Property = this.props.getProperty(PROPS.rowTitleWidth);
        rowTitleWidthProp.options.max = cw;
      });
    });
  }

  // -- Routines --


  // -- Handle events --

  private rootOnSizeChanged(_value: any): void {
    _package.logger.debug('BiotrackViewer.rootOnSizeChanged() ');
    this.calcSize();
  }

  // -- IRenderer--

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    // Put the event trigger in the tail of the synced calls queue.
    this.viewSyncer.sync('invalidate(${caller ? ` <- ${caller} ` : \'\'})', async () => {
      // update view / render
      this._onRendered.next();
    });
  }

  async awaitRendered(timeout: number | undefined = 5000): Promise<void> {
    await testEvent(this.onRendered, () => {}, () => {
      this.invalidate();
    }, timeout);

    // Rethrow stored syncer error (for test purposes)
    const viewErrors = this.viewSyncer.resetErrors();
    if (viewErrors.length > 0) throw viewErrors[0];
  }
}
