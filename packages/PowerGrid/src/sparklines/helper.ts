// need for calculate distance between mouse and point
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export function distance(p1: DG.Point, p2: DG.Point): number {
  return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

export class Hit {
  activeColumn: number = -1;
  cols: DG.Column[] = [];
  row: number = -1;
  isHit: boolean = false;
}

export function createTooltip(cols: any, activeColumn: number, row: any): any {
  let arr = [];
  for (let i = 0; i < cols.length; i++) {
    arr.push(ui.divH([ui.divText(`${cols[i].name}:`, {
          style: {
            margin: '0 10px 0 0',
            fontWeight: (activeColumn == i) ? 'bold' : 'normal',
          }
        }), ui.divText(`${Math.floor(cols[i].get(row) * 100) / 100}`, {
          style: {
            fontWeight: (activeColumn == i) ? 'bold' : 'normal',
          }
        })]
      )
    );
  }
  return arr;
}
