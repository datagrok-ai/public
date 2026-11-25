import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {syncItemWithAxes} from './formula-lines/axis-sync';
import {CreationControl} from './formula-lines/creation-control';
import {Editor} from './formula-lines/editor';
import {Host} from './formula-lines/host';
import {ItemList} from './formula-lines/item-list';
import {Preview, supportsRegionDrawing} from './formula-lines/preview';
import {Table} from './formula-lines/table';
import {DEFAULT_OPTIONS, EditorItem, EditorOptions, ITEM_SOURCE} from './formula-lines/types';

export {DEFAULT_OPTIONS};
export type {EditorOptions};

/**
 * A Dialog window with Formula Lines list, preview and editor.
 */
export class FormulaLinesDialog {
  public dialog: DG.Dialog = ui.dialog({
    title: 'Formula Lines',
    helpUrl: '/help/develop/how-to/viewers/show-formula-lines.md',
  });

  private host: Host;
  private preview: Preview;
  private editor: Editor;
  private viewerTable?: Table;
  private dframeTable?: Table;
  private creationControl: CreationControl;
  private tabs: DG.TabControl;
  // Items that failed the validation of their form; OK is disabled while there are any
  private invalidItems = new Set<EditorItem>();

  /** Returns the Table corresponding to the current tab in the tab control */
  public get currentTable(): Table {
    return this.tabs.currentPane.name === ITEM_SOURCE.VIEWER ? this.viewerTable! : this.dframeTable!;
  }

  /** Initializes all parameters and opens the Dialog window */
  constructor(
    src: DG.DataFrame | DG.Viewer,
    private options: EditorOptions = DEFAULT_OPTIONS,
    private showValueOnOpen?: { index?: number, isDataFrame?: boolean, isAnnotationArea?: boolean }) {
    /** Init Helpers */
    this.host = new Host(src);
    this.creationControl = this.initCreationControl();
    this.preview = this.initPreview();
    this.editor = this.initEditor();
    this.tabs = this.initTabs();
    this.dialog.sub(this.dialog.onClose.subscribe(() => this.dialog.detach()));

    const width = Math.min(1000, Math.floor(document.body.clientWidth / 1.3));
    const height = Math.min(800, Math.floor(document.body.clientHeight / 1.5));
    this.dialog
      .add(this.createLayout())
      .onOK(this.onOKAction.bind(this), {closeOnEnter: false})
      .show({
        resizable: true,
        width,
        height,
        x: Math.floor((window.innerWidth - width) / 2),
        y: Math.floor((window.innerHeight - height) / 2),
      });

    this.initDefaultOnOpenState();

    this.dialog.sub(syncItemWithAxes(this.preview, this.editor, () => this.currentTable));
  }

  // The items table over the preview on the left, the editor on the right
  private createLayout(): HTMLElement {
    const tabsInitialHeight = '140px';
    const editorInitialWidth = '420px';
    const makeWrap = (child: HTMLElement): HTMLDivElement => {
      const wrapDiv = document.createElement('div');
      // plain <div>, not ui.div — ui.div tags it with `ui-div`, which would trigger
      // `div.ui-div > div.ui-box { width:400px; height:300px }` on the d4-viewer/d4-tab-host child
      wrapDiv.className = 'ui-update-shadow dlg-formula-lines-pane-wrap';
      wrapDiv.appendChild(child);
      return wrapDiv;
    };

    const tabsBox = ui.box(makeWrap(this.tabs.root),
      {classes: 'dlg-formula-lines-tabs-box', style: {height: tabsInitialHeight}});

    const previewBox = ui.box(makeWrap(this.preview.root), {classes: 'dlg-formula-lines-preview-box'});
    const leftSplit = ui.splitV([tabsBox, previewBox], {classes: 'dlg-formula-lines-left-split'}, true);
    const editorBox = ui.box(this.editor.root, {classes: 'dlg-formula-lines-editor-box',
      style: {width: editorInitialWidth}});

    return ui.splitH([leftSplit, editorBox],
      {style: {width: '100%', height: '100%', minHeight: '0', minWidth: '0'}}, true);
  }

  private initDefaultOnOpenState(): void {
    if (!this.showValueOnOpen)
      return;

    // A table has no Viewer tab, and a viewer may come without the DataFrame one: the first tab stays then
    const pane = this.tabs.getPane(this.showValueOnOpen.isDataFrame ? ITEM_SOURCE.DATAFRAME : ITEM_SOURCE.VIEWER);
    if (pane)
      this.tabs.currentPane = pane;
    const table = this.currentTable;
    if (this.showValueOnOpen.index != null)
      table.update(table.items.rowIdx(this.showValueOnOpen.index, !this.showValueOnOpen.isAnnotationArea));
    else
      table.setFirstItemAsCurrent();
  }

  private initPreview(): Preview {
    const preview = new Preview(this.host.initialItems, this.host.src, this.creationControl.popupMenu);
    preview.otherFormulaLines = this.host.otherFormulaLines(this.host.initialItems);
    return preview;
  }

  private initCreationControl(): CreationControl {
    return new CreationControl({
      getCols: () => this.preview.axisCols,
      getCurrentItem: () => this.currentTable.currentItem,
      onItemCreated: (item) => this.onItemCreatedAction(item),
      createArea: (lassoMode) => this.drawArea(lassoMode),
    }, this.host);
  }

  // Resolves to null if the preview viewer has no drawing
  private drawArea(lassoMode?: boolean): Promise<DG.AnnotationRegion | null> {
    const viewer = this.preview.viewer;
    if (!supportsRegionDrawing(viewer))
      return Promise.resolve(null);

    return new Promise((resolve) => {
      viewer.disableAnnotationRegionDrawing();
      viewer.setOptions({annotationRegions: '[]', formulaLines: '[]'});
      this.editor.update(-1);
      viewer.enableAnnotationRegionDrawing(lassoMode, (region: { [key: string]: unknown }) => {
        region['isDataFrameRegion'] = this.tabs.currentPane.name === ITEM_SOURCE.DATAFRAME;
        const props = viewer.props as DG.IScatterPlotSettings;
        const annotationRegions = JSON.parse(props.annotationRegions || '[]');
        annotationRegions.push(region);
        props.annotationRegions = JSON.stringify(annotationRegions);
        resolve(region as DG.AnnotationRegion);
      });
    });
  }

  private initEditor(): Editor {
    return new Editor(this.host.initialItems, this.preview.dataFrame,
      (rowIdx: number): boolean => {
        this.currentTable.update(rowIdx);
        return this.preview.update(rowIdx);
      },
      (item, isValid) => this.setItemValidity(item, isValid),
      this.host.singleAxisColumnName,
      this.host.defaultItemZIndex,
    );
  }

  private setItemValidity(item: EditorItem, isValid: boolean): void {
    isValid ? this.invalidItems.delete(item) : this.invalidItems.add(item);
    const ok = this.dialog.getButton('OK');
    ok.classList.toggle('disabled', this.invalidItems.size > 0);
    this.invalidItems.size > 0 ? ok.setAttribute('aria-disabled', 'true') : ok.removeAttribute('aria-disabled');
  }

  private initTabs(): DG.TabControl {
    const tabs = DG.TabControl.create();

    /** Init Viewer Table (in the first tab) */
    if (this.host.viewerItems) {
      tabs.addPane(ITEM_SOURCE.VIEWER, () => {
        this.viewerTable = this.initTable(this.host.viewerItems!);
        return this.viewerTable.root;
      });
    }

    /** Init DataFrame Table (in the second tab) */
    if (this.options.allowEditDFLines) {
      tabs.addPane(ITEM_SOURCE.DATAFRAME, () => {
        this.dframeTable = this.initTable(this.host.dframeItems);
        return this.dframeTable.root;
      });
    }
    // Overrides the standard component logic that hides the header containing only one tab
    tabs.header.style.removeProperty('display');

    /** Display "Add new" button */
    tabs.header.append(this.creationControl.button);

    /** Change data source when switching tabs */
    tabs.onTabChanged.subscribe((_) => {
      const table = this.currentTable;
      this.editor.items = table.items;
      this.preview.items = table.items;
      this.preview.otherFormulaLines = this.host.otherFormulaLines(table.items);
      table.setFirstItemAsCurrent();
    });

    return tabs;
  }

  private initTable(items: ItemList): Table {
    return new Table(items,
      (rowIdx: number): void => {
        this.editor.update(rowIdx);
        this.preview.update(rowIdx);
      }, (item) => this.setItemValidity(item, true), this.preview.srcAxes, !this.showValueOnOpen);
  }

  private onOKAction() {
    this.host.save();
    this.creationControl.saveHistory();
  }

  private onItemCreatedAction(item: EditorItem): void {
    const rowIdx = this.currentTable.add(item);
    this.editor.update(rowIdx);
    this.preview.update(rowIdx);
  }
}
