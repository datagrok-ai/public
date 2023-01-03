import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {MarkupNodeType, TreeStylerBase} from './markup';
import {isLeaf} from '@datagrok-libraries/bio';


export class DendrogramTreeStyler extends TreeStylerBase<MarkupNodeType> {
  override get lineWidth(): number { return this._lineWidth; }

  set lineWidth(value: number) {
    this._lineWidth = value;
    this._onStylingChanged.next();
  }

  override get nodeSize(): number { return this._nodeSize; }

  set nodeSize(value: number) {
    this._nodeSize = value;
    this._onStylingChanged.next();
  }

  override get showGrid(): boolean { return this._showGrid; }

  set showGrid(value: boolean) {
    this._showGrid = value;
    this._onStylingChanged.next();
  }

  get strokeColor(): string { return this._strokeColor; }

  set strokeColor(value: string) {
    this._strokeColor = value;
    this._onStylingChanged.next();
  }

  override getStrokeColor(node: MarkupNodeType): string { return this._strokeColor; }

  get fillColor(): string { return this._fillColor; }

  set fillColor(value: string) {
    this._fillColor = value;
    this._onStylingChanged.next();
  }

  override getFillColor(node: MarkupNodeType): string { return this._fillColor; }

  constructor(name: string,
    lineWidth: number, nodeSize: number, showGrid: boolean, strokeColor?: string, fillColor?: string
  ) {
    super(name, lineWidth, nodeSize, showGrid, strokeColor, fillColor);
  }
}

export class DendrogramColorCodingTreeStyler extends DendrogramTreeStyler {
  protected _nodeCol: DG.Column;

  protected _colorCol: DG.Column;

  protected _colorAggrType: string;

  protected _rowByNameDict: { [nodeName: string]: number };

  constructor(name: string, lineWidth: number, nodeSize: number, showGrid: boolean,
    nodeCol: DG.Column, colorCol: DG.Column, colorAggrType: string,
    strokeColor: string, fillColor: string
  ) {
    super(name, lineWidth, nodeSize, showGrid, strokeColor, fillColor);
    this._nodeCol = nodeCol;
    this._colorCol = colorCol;
    this._colorAggrType = colorAggrType;

    this._rowByNameDict = {};
    const dfLength = this._colorCol.length;
    for (let rowI: number = 0; rowI < dfLength; rowI++) {
      const nodeName: string = this._nodeCol.get(rowI);
      this._rowByNameDict[nodeName] = rowI;
    }
  }

  override getStrokeColor(node: MarkupNodeType): string {
    let res: string;
    const colorType: string = this._colorCol.colors.getType();
    if (colorType != DG.COLOR_CODING_TYPE.OFF) {
      const nodeName: string = node.name;
      const rowI: number = this._rowByNameDict[nodeName];
      const color: number = this._colorCol.colors.getColor(rowI);
      res = DG.Color.toRgb(color);
    } else {
      res = this._strokeColor;
    }
    return res;
  }

  override getFillColor(node: MarkupNodeType): string {
    let res: string;
    const colorType: string = this._colorCol.colors.getType();
    if (colorType != DG.COLOR_CODING_TYPE.OFF) {
      // const val = this._colorCol.get(node.index);
      const nodeName: string = node.name;
      const rowI: number = this._rowByNameDict[nodeName];
      const color: number = this._colorCol.colors.getColor(rowI);
      res = DG.Color.toRgb(color);
    } else {
      res = this._fillColor;
    }
    return res;
  }
}
