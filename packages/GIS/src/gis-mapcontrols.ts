//Base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {OpenLayers} from '../src/gis-openlayer';
import {Control, defaults as defaultControls} from 'ol/control'; //Attribution include?


export class PanelLayersControl extends Control {
  /**
   * @param {Object} [opt] Control options.
   */
  // parentOnClickHandler: Function | null = null;
  dfLayersList: DG.DataFrame | undefined;
  layersGrid: DG.Grid;
  layersPanel: HTMLElement;
  ol: OpenLayers;

  constructor(parent: OpenLayers, opt: any | undefined, clkHandler: Function | null = null) {
    const options = opt || {};

    const layersGridStyle = {
      'allowEdit': false,
      'showAddNewRowIcon': true,
      'allowColReordering': false,
      'allowColResizing': false,
      'allowColHeaderResizing': false,
      'allowColSelection': false,
      'allowRowResizing': false,
      'allowRowSelection': false,
      'allowBlockSelection': false,
      'autoScrollColumnIntoView': false,
      'showCurrentRowIndicator': false,
      'showCurrentCellOutline': false,
      'showRowHeader': false,
      'showColumnTooltip': false,
      'showColumnGridlines': false,
      'topLevelDefaultMenu': true,
      'showColumnLabels': false,
    };

    const df = grok.data.demo.demog(6);
    // const layersGrid = DG.Viewer.grid(this.dfLayersList, layersGridStyle);
    const layersGrid = DG.Viewer.grid(df, layersGridStyle);
    layersGrid.autoSize(200, 170);

    // const element = ui.box(layersGrid.root);
    const element = ui.div(layersGrid.root);
    // const element = ui.div();
    element.className = 'panel-layers ol-unselectable ol-control';

    super({
      element: element,
      target: options.target,
    });

    this.ol = parent;
    this.layersGrid = layersGrid;
    this.layersPanel = element;
    this.refreshDF(DG.DataFrame.fromColumns([]));
    //<<PanelLayersControl constructor
  }

  setupGridControls() {
    this.layersGrid.setOptions({'rowHeight': 24});
    // this.layersGrid.rows
    let col = this.layersGrid.columns.byName('vis');
    if (col) {
      col.width = 19;
      col.cellType = 'html';
      col.editable = true;
    }
    col = this.layersGrid.columns.byName('del');
    if (col) {
      col.width = 19;
      col.cellType = 'html';
    }
    col = this.layersGrid.columns.byName('exp');
    if (col) {
      col.width = 19;
      col.cellType = 'html';
    }
    col = this.layersGrid.columns.byName('layerid');
    if (col)
      col.visible = false;
    col = this.layersGrid.columns.byName('name');
    if (col)
      col.width = 125;

    this.layersGrid.onCellPrepare(function(gc) {
      if (gc.isTableCell && gc.gridColumn.name === 'vis') {
        const btnVisible = ui.iconFA(gc.cell.value ? 'eye' : 'eye-slash', undefined, 'Show/hide layer');

/*        const btnVisible = ui.iconFA(gc.cell.value ? 'eye' : 'eye-slash', (evt) => {
          // evt.stopImmediatePropagation();
          //get layer ID for grid row >>
          const layerId = 12;
          if (!layerId)
            return;
          // const layer = this.ol.getLayerById(layerId);
          // if (!layer)
          //   return;
          // const isVisible = layer.getVisible();
          // layer.setVisible(!isVisible);
          // this.updateLayersList();
        }, 'Show/hide layer'); */
        btnVisible.setAttribute('LayerId', 'testID');
        btnVisible.style.margin = '3px';

        gc.style.element = ui.divV([
          btnVisible,
        ]);
      } //setup visible button

      if (gc.isTableCell && gc.gridColumn.name === 'del') {
        const btnDel = ui.iconFA('trash-alt', (evt) => {
        }, 'Delete layer');
        btnDel.style.margin = '3px';
        gc.style.element = ui.divV([btnDel]);
      } //setup delete button

      if (gc.isTableCell && gc.gridColumn.name === 'exp') {
        const btnExp = ui.iconFA(gc.cell.value ? 'download' : '', (evt) => {
        }, 'Export layer data to table');
        btnExp.style.margin = '3px';
        gc.style.element = ui.divV([btnExp]);
      } //setup export button
    });
  }

  updateLayersList() {
    if (this.ol)
      this.ol.updateLayersList();
  }

  refreshDF(df: DG.DataFrame | undefined) {
    if (!df)
      return;
    this.dfLayersList = df;
    this.layersGrid.dataFrame = this.dfLayersList;
    this.setupGridControls();
    this.layersGrid.invalidate();
  }
  setVisibility(visible: boolean) {
    // this.setProperties({'visible': visible});
    this.layersPanel.style.visibility = visible ? 'hidden' : 'visible';
  }
}

export class BtnLayersControl extends Control {
  /**
   * @param {Object} [opt] Control options.
   */
  parentOnClickHandler: Function | null = null;
  btnIsOn = false;

  constructor(opt: any | undefined, clkHandler: Function | null = null) {
    const options = opt || {};

    // const button = document.createElement('button');
    const button = ui.iconFA('layer-group');
    // button.innerHTML = 'L';
    // button.className = 'btn-layers ol-unselectable ol-control';

    const element = document.createElement('div');
    // const element = ui.iconFA('layer-group');
    element.className = 'btn-layers ol-unselectable ol-control';
    element.appendChild(button);

    super({
      element: element,
      target: options.target,
    });

    this.parentOnClickHandler = clkHandler;
    button.addEventListener('click', this.onClickLayersBtnControl.bind(this), false);
  }

  onClickLayersBtnControl() {
    this.btnIsOn = !this.btnIsOn;
    if (this.parentOnClickHandler)
      this.parentOnClickHandler();
  }
}