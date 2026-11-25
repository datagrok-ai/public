import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {formatAreaFormula} from './formulas';
import {ItemHistory} from './history';
import {createItem, NewItemDefaults} from './item-factory';
import {AxisColumns, BTN_CAPTION, EditorItem, isAnnotationRegionType, ITEM_CAPTION, ITEM_TYPE} from './types';

const HISTORY_KEY = 'formula-lines-dialog';
const HISTORY_KEY_ANNOTATIONS = 'annotation-regions-dialog';

/** What the creation of an item needs from the dialog. */
export interface CreationActions {
  /** The columns on the axes of the preview: the new items refer to them */
  getCols: () => AxisColumns,
  /** The item to clone */
  getCurrentItem: () => EditorItem | null,
  /** Updates the Table, Preview and Editor states after item creation */
  onItemCreated: (item: EditorItem) => void,
  /** Lets the user draw an area on the preview; null if the preview viewer has no drawing */
  createArea: (lassoMode?: boolean) => Promise<DG.AnnotationRegion | null>,
}

/** What the host viewer allows to create, and its defaults for the new items. */
export interface CreationOptions extends NewItemDefaults {
  /** Allow-list of `ITEM_CAPTION` strings to show in the popup. `undefined` means "all five". */
  allowedItems?: string[],
}

/**
 * Helper that implements the logic of creating a Formula Line item of a given type.
 */
export class CreationControl {
  private formulaLinesHistory = new ItemHistory<DG.FormulaLine>(HISTORY_KEY);
  private annotationRegionsHistory = new ItemHistory<DG.AnnotationRegion>(HISTORY_KEY_ANNOTATIONS);

  constructor(private actions: CreationActions, private options: CreationOptions) {}

  /** Creates a button and binds an item creation menu to it */
  public get button(): HTMLElement {
    const btn = ui.bigButton(BTN_CAPTION.ADD_NEW, this.popupMenu);
    return ui.div([btn], {style: {width: '100%', textAlign: 'right'}});
  }

  public saveHistory() {
    this.formulaLinesHistory.save();
    this.annotationRegionsHistory.save();
  }

  /** Opens a popup menu with predefined new Formula Line item types.
   *  [valY] and [valX] are the point of the preview the menu is opened at: constant lines pass through it. */
  public popupMenu = (valY?: number, valX?: number): void => {
    const onClickAction = (itemCaption: string) => this.createItem(itemCaption, valY, valX);
    const {allowedItems} = this.options;

    /** Construct popup menu — filtered by host viewer's allow-list when set */
    const allItems = [
      ITEM_CAPTION.LINE,
      ITEM_CAPTION.VERT_LINE,
      ITEM_CAPTION.HORZ_LINE,
      ITEM_CAPTION.VERT_BAND,
      ITEM_CAPTION.HORZ_BAND,
    ];
    const items = allowedItems ? allItems.filter((c) => allowedItems.includes(c)) : allItems;
    const menu = DG.Menu.popup().items(items, onClickAction);

    this.addRegionItems(menu, onClickAction);

    const currentItem = this.actions.getCurrentItem();
    const formulaLines = this.formulaLinesHistory.items;
    const annotationRegions = this.annotationRegionsHistory.items;

    /** Add separator only if other menu items exist */
    if (currentItem || formulaLines.length || annotationRegions.length)
      menu.separator();

    /**
     * Add "Clone" menu if the current table line exists.
     * TODO: The best option is to make the menu item enabled/disabled. But there is no such API yet.
     */
    if (currentItem)
      menu.items([BTN_CAPTION.CLONE], onClickAction);

    /**
     * Add "History" menu group.
     * TODO: The best option is to make the menu item enabled/disabled. But there is no such API yet.
     */
    if (formulaLines.length) {
      this.fillHistoryGroup(menu.group(BTN_CAPTION.FORMULA_LINES_HISTORY), this.formulaLinesHistory,
        (item) => item.formula!);
    }

    if (annotationRegions.length) {
      this.fillHistoryGroup(menu.group(BTN_CAPTION.ANNOTATION_REGIONS_HISTORY), this.annotationRegionsHistory,
        formatAreaFormula);
    }

    menu.show();
  };

  // The region items the axes of the preview can show
  private addRegionItems(menu: DG.Menu, onClickAction: (itemCaption: string) => void): void {
    const {allowedItems} = this.options;
    const addItem = (itemCaption: string, description: string) =>
      menu.item(itemCaption, () => onClickAction(itemCaption), null, {description});

    const cols = this.actions.getCols();
    // Bands don't use axis maps; any xMap/yMap leaking through from the preview viewer is spurious.
    const currentItem = this.actions.getCurrentItem();
    const isBandCurrent = currentItem && !isAnnotationRegionType(currentItem.type) &&
      (currentItem as DG.FormulaLine).type === ITEM_TYPE.BAND;
    const xMapActive = !isBandCurrent && !!cols.xMap;
    const yMapActive = !isBandCurrent && !!cols.yMap;
    // 2D hosts: full set of region items (formula, draw rect, draw lasso).
    // 1D hosts (histogram, bar chart, box plot — `allowedItems != null`) get a draw-rectangle
    // option that the Dart-side `createDrawnRegion` override converts into a value-axis formula
    // region, plus a single axis-bounds range item targeting the value column.
    const xRegionEligible = !!cols.x?.isNumerical && !xMapActive;
    const yRegionEligible = !!cols.y?.isNumerical && !yMapActive;
    if (allowedItems == null) {
      if (xRegionEligible && yRegionEligible) {
        addItem(ITEM_CAPTION.FORMULA_REGION, 'Adds a new area defined by two formula lines.');
        addItem(ITEM_CAPTION.RECT_REGION, 'Draws a rectangle area.');
        addItem(ITEM_CAPTION.POLYGON_REGION, 'Draws a polygon area using lasso tool.');
      }
      return;
    }

    // Single-axis host. Value lives on the X axis when `allowedItems` includes vertical
    // bands (histogram, horizontal bar chart) and on the Y axis otherwise (vertical bar
    // chart, box plot).
    const valueOnX = allowedItems.includes(ITEM_CAPTION.VERT_BAND);
    const rangeDescription = 'Adds a region bounded by two values on the value axis (q1 .. q3 of the column).';
    if (xRegionEligible || yRegionEligible)
      addItem(ITEM_CAPTION.RECT_REGION, 'Draws a rectangle and converts it to a value-axis formula region.');
    if (valueOnX && xRegionEligible)
      addItem(ITEM_CAPTION.VERT_RANGE_REGION, rangeDescription);
    else if (!valueOnX && yRegionEligible)
      addItem(ITEM_CAPTION.HORZ_RANGE_REGION, rangeDescription);
  }

  private fillHistoryGroup<T extends EditorItem>(group: DG.Menu, history: ItemHistory<T>,
    getTitle: (item: T) => string): void {
    for (const item of history.items) {
      const title = getTitle(item);
      if (title) {
        group.item(title.length > 50 ? title.substring(0, 50) + '...' : title, () => {
          const newItem = structuredClone(item);
          history.add(newItem);
          this.actions.onItemCreated(newItem);
        });
      }
    }

    group.endGroup();
  }

  private createItem(itemCaption: string, valY?: number, valX?: number): void {
    if (itemCaption === ITEM_CAPTION.POLYGON_REGION || itemCaption === ITEM_CAPTION.RECT_REGION) {
      this.actions.createArea(itemCaption === ITEM_CAPTION.POLYGON_REGION).then((region) => {
        if (!region)
          return;

        region.zIndex ??= this.options.defaultItemZIndex;
        this.onCreated(region);
      });
      return;
    }

    const currentItem = itemCaption === BTN_CAPTION.CLONE ? this.actions.getCurrentItem() : null;
    const item = createItem(itemCaption, this.actions.getCols(), this.options, currentItem, valY, valX);
    if (item)
      this.onCreated(item);
  }

  private onCreated(item: EditorItem): void {
    (isAnnotationRegionType(item.type) ? this.annotationRegionsHistory : this.formulaLinesHistory).add(item);
    this.actions.onItemCreated(item);
  }
}
