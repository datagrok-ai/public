import * as DG from 'datagrok-api/dg';


const minLogFloat = 1e-30;

/** Represents 1D transformation between screen and world coordinates */
interface ITransform {
  screenToWorld(screen: number, length: number, min: number, max: number, inverse?: boolean): number;
  worldToScreen(world: number, length: number, min: number, max: number, inverse?: boolean): number;
}


/** Linear transformation */
class _LinearTransform implements ITransform {
  screenToWorld(screen: number, length: number, min: number, max: number, inverse?: boolean): number {
    const d = inverse ? 1 - screen / length: screen / length;
    return min + d * (max - min);
  }

  worldToScreen(world: number, length: number, min: number, max: number, inverse?: boolean): number {
    return inverse ?
      ((1 - ((world - min)) / (max - min)) * length) :
      ((world - min) / (max - min) * length);
  }
}


/** Log transformation */
class _LogTransform implements ITransform {
  screenToWorld(screen: number, length: number, min: number, max: number, inverse?: boolean): number {
    min = (min < minLogFloat ? minLogFloat : min);

    const maxLog = Math.log(max);
    const minLog = Math.log(min);
    return Math.exp(minLog + (maxLog - minLog) * (inverse ? length - screen : screen) / length);
  }

  worldToScreen(world: number, length: number, min: number, max: number, inverse?: boolean): number {
    if (world <= 0) return NaN;

    //eliminates problems with rounding
    if (Math.abs(world - min) < minLogFloat) return inverse ? length : 0.0;
    if (Math.abs(world - max) < minLogFloat) return inverse ? 0.0 : length;

    min = (min < minLogFloat ? minLogFloat : min);
    world = (world < minLogFloat ? minLogFloat : world);

    const res = (length * Math.log(world / min) / Math.log(max / min));
    return inverse ? length - res : res;
  }
}

const linearTransform = new _LinearTransform();
const logTransform = new _LogTransform();


/** Provides convenient methods for mapping between screen and world coordinates */
export class Viewport {
  world: DG.Rect;
  screen: DG.Rect;
  logX: boolean = false;
  logY: boolean = false;
  inverseX: boolean = false;
  inverseY: boolean = false;

  get xt(): ITransform { return this.logX ? logTransform : linearTransform; }
  get yt(): ITransform { return this.logY ? logTransform : linearTransform; }

  constructor(world: DG.Rect, screen: DG.Rect, logX = false, logY = false) {
    this.world = world;
    this.screen = screen;
    this.logX = logX;
    this.logY = logY;
  }

  xToScreen(world: number): number {
    return this.screen.left + this.xt.worldToScreen(world, this.screen.width,
      this.world.minX, this.world.maxX, this.inverseX);
  }

  yToScreen(world: number): number {
    return this.screen.top + this.yt.worldToScreen(world, this.screen.height,
      this.world.minY, this.world.maxY, !this.inverseY);
  }

  xToWorld(screen: number): number {
    return this.xt.screenToWorld(screen - this.screen.left, this.screen.width,
      this.world.minX, this.world.maxX, this.inverseX);
  }

  yToWorld(screen: number): number {
    return this.yt.screenToWorld(screen - this.screen.top, this.screen.height,
      this.world.minY, this.world.maxY, !this.inverseY);
  }

  toScreen(world: DG.Point): DG.Point { return new DG.Point(this.xToScreen(world.x), this.yToScreen(world.y)); }
  toWorld(screen: DG.Point): DG.Point { return new DG.Point(this.xToWorld(screen.x), this.yToWorld(screen.y)); }

  drawCoordinateGrid(g: CanvasRenderingContext2D, xAxis: DG.Rect | undefined, yAxis: DG.Rect | undefined): void {
    if (xAxis) {
      DG.Paint.horzAxis(g, this.world.minX, this.world.maxX, xAxis.x, xAxis.y,
        xAxis.width, xAxis.height, this.logX, this.inverseX);
    }
    if (yAxis) {
      DG.Paint.vertAxis(g, this.world.minY, this.world.maxY, yAxis.x, yAxis.y,
        yAxis.width, yAxis.height, this.logY, this.inverseY);
    }
  }
}
