//Base import
//import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//import $ from 'cash-dom';
import {OpenLayers} from '../src/gis-openlayer';
import {Control} from 'ol/control'; //Attribution include?


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
      // 'allowEdit': false,
      'showAddNewRowIcon': false,
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
      'showColumnLabels': false,
      'topLevelDefaultMenu': true,
      'allowDynamicMenus': false,
      'showContextMenu': false,
    };

    const df = DG.DataFrame.fromCsv(
      `vis, name, exp, del, layerid
      true, Map, false, false, 0
      `);
    const layersGrid = DG.Viewer.grid(df, layersGridStyle);
    layersGrid.autoSize(200, 300, 200, 100, true);
    layersGrid.columns.setOrder(['vis', 'name', 'exp', 'del']);
    layersGrid.root.style.borderWidth = '2px';
    layersGrid.root.style.borderColor = 'rgb(100, 100, 100)';

    const element = ui.div(layersGrid.root);
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

  cellPrepareFn(gc: DG.GridCell) {
    if (gc.isTableCell && gc.gridColumn.name === 'del') {
      const btnDel = ui.button(ui.iconFA(gc.cell.value ? 'trash-alt' : ''), () => {
        if (gc.tableRowIndex !== null) {
          const col = gc.grid.table.columns.byName('layerid');
          if (col) {
            const layerId = col.get(gc.tableRowIndex);
            if (layerId)
              this.ol.removeLayerById(layerId);
            this.updateLayersList();
          }
        }
      }, 'Delete layer');
      btnDel.style.width = '18px';
      btnDel.style.height = '18px';
      gc.style.element = btnDel;
    } //setup delete button

    if (gc.isTableCell && gc.gridColumn.name === 'exp') {
      const btnExp = ui.button(ui.iconFA(gc.cell.value ? 'arrow-to-bottom' : ''), () => {
        if (gc.tableRowIndex !== null) {
          const col = gc.grid.table.columns.byName('layerid');
          if (col) {
            const layerId = col.get(gc.tableRowIndex);
            // if (layerId)
            //   this.ol.removeLayerById(layerId);
            this.updateLayersList();
          }
        }
      }, 'Export layer data to table');
      btnExp.style.width = '18px';
      btnExp.style.height = '18px';
      gc.style.element = btnExp;
    } //setup export button
  }

  cellVisibleClick(gc: DG.GridCell) {
    if (gc.tableRowIndex !== null) {
      const col = gc.grid.table.columns.byName('layerid');
      if (col) {
        const layerId = col.get(gc.tableRowIndex);
        if (!layerId)
          return;
        const layer = this.ol.getLayerById(layerId);
        if (!layer)
          return;

        layer.setVisible(gc.cell.value);
        this.updateLayersList();
      }
    }
  }

  setupGridControls() {
    this.layersGrid.setOptions({'rowHeight': 24});
    let col = this.layersGrid.columns.byName('vis');
    if (col) {
      col.width = 19;
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
    if (col) {
      col.width = 125;
      col.editable = false;
    }

    this.layersGrid.onCellPrepare(this.cellPrepareFn.bind(this));
    this.layersGrid.onCurrentCellChanged.subscribe(this.cellVisibleClick.bind(this));
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
    if (visible)
      this.updateLayersList();
    this.layersPanel.style.visibility = visible ? 'visible' : 'hidden';
  }
}

export class BtnLayersControl extends Control {
  /**
   * @param {Object} [opt] Control options.
   */
  parentOnClickHandler: Function | null = null;
  isOn = false;

  constructor(opt: any | undefined, clkHandler: Function | null = null) {
    const options = opt || {};

    const button = ui.button(ui.iconFA('layer-group'), ()=>{}, 'Map layers');
    button.className = 'btn-layers ol-unselectable ol-control';

    const element = document.createElement('div');
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
    this.isOn = !this.isOn;
    if (this.parentOnClickHandler)
      this.parentOnClickHandler();
  }
}

//lPanel: PanelLayersControl;
