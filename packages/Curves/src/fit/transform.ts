import * as DG from 'datagrok-api/dg';

export interface ITransform {
  xToScreen(world: number): number;
  yToScreen(world: number): number;
}


export class Transform {
  static linear(world: DG.Rect, screen: DG.Rect): ITransform {
    return {
      xToScreen(worldX: number): number {
        return screen.left + screen.width * (worldX - world.left) / world.width;
      },
      yToScreen(worldY: number): number {
        return screen.bottom - screen.height * (worldY - world.top) / world.height;
      }
    };
  }
}


export class CanvasChart {
  g: CanvasRenderingContext2D;
  world: DG.Rect = new DG.Rect(0, 0, 1, 1);
  screen: DG.Rect;
  logX: boolean = false;
  logY: boolean = false;

  xToScreen(world: number): number { return 0; }
  yToScreen(world: number): number { return 0; }

  constructor(g: CanvasRenderingContext2D, screen: DG.Rect) {
    this.g = g;
    this.screen = screen;
  }
}