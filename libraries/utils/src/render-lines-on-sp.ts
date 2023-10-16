import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import BitArray from './bit-array';

export type ILines = {
    from: Uint32Array;
    to: Uint32Array;
    colors?: string[];
    widths?: Uint32Array;
    opacities?: Float32Array;
    drawArrows?: BitArray;
}

export class ScatterPlotLinesRenderer {
    sp: DG.ScatterPlotViewer;
    xAxisValues: any[];
    yAxisValues: any[];
    currentLineIdx = -1;
    lines: ILines;
    lineClicked = new Subject<number>();
    lineHover = new Subject<number>();
    canvas: HTMLCanvasElement;
    ctx: CanvasRenderingContext2D;
    paths: Path2D[];
    mouseOverLineId: number | null = null;
    numberOfLinesPerPair: {[key: string]: number} = {};

    get currentLineId(): number {
        return this.currentLineIdx;
    }

    set currentLineId(id: number) {
        this.currentLineIdx = id;
        this.sp.render(this.canvas.getContext('2d') as CanvasRenderingContext2D);
    }

    constructor(sp: DG.ScatterPlotViewer, xAxis: string, yAxis: string, lines: ILines) {
        this.sp = sp;
        this.xAxisValues = this.sp.dataFrame!.columns.byName(xAxis).toList();
        this.yAxisValues = this.sp.dataFrame!.columns.byName(yAxis).toList();;
        this.lines = lines;
        this.canvas = this.sp.getInfo()['canvas'];
        this.ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        this.paths = new Array<Path2D>(this.lines.from.length);
        this.countLinesPerPair();

        this.canvas.onmousedown = () => {
            if (this.mouseOverLineId !== null)
                this.lineClicked.next(this.mouseOverLineId);
        }

        this.canvas.onmousemove = (event: MouseEvent) => {
            this.mouseOverLineId = this.checkCoordsOnLine(event.offsetX, event.offsetY);
            if (this.mouseOverLineId !== null)
                this.lineHover.next(this.mouseOverLineId);
        }

        sp.onEvent('d4-before-draw-scene')
            .subscribe((_: any) => {
                this.renderLines();
            });
    }


    renderLines(): void {
        const filterBitset = this.sp.dataFrame.filter;
        const linesPerPairCounter: { [key: string]: number } = {};
        for (let i = 0; i < this.lines.from.length; i++) {
            if (filterBitset.get(this.lines.from[i]) && filterBitset.get(this.lines.to[i])) {
                const pointFrom = this.sp.worldToScreen(this.xAxisValues[this.lines.from[i]], this.yAxisValues[this.lines.from[i]]);
                const aX = pointFrom?.x;
                const aY = pointFrom?.y;
                const pointTo = this.sp.worldToScreen(this.xAxisValues[this.lines.to[i]], this.yAxisValues[this.lines.to[i]]);
                const bX = pointTo?.x;
                const bY = pointTo?.y;
                if (aX && aY && bX && bY) {
                    const pairStr = this.getStringPair(this.lines.from[i], this.lines.to[i]);
                    !linesPerPairCounter[pairStr] ? linesPerPairCounter[pairStr] = 1 : linesPerPairCounter[pairStr]++;
                    const line = new Path2D();
                    this.paths[i] = line;
                    line.moveTo(aX!, aY!);
                    const color = this.lines.colors?.[i] ? this.lines.colors?.[i] : '0,128,0';
                    const opacity = this.lines.opacities?.[i] ? this.lines.opacities?.[i] : 1;
                    this.ctx.strokeStyle = `rgba(${color},${opacity})`;
                    this.ctx.lineWidth = this.lines.widths?.[i] ? i === this.currentLineIdx ? this.lines.widths?.[i] + 2 :
                    this.lines.widths?.[i] : i === this.currentLineIdx ? 3 : 1;
                    let midPoint: DG.Point = this.midPoint(pointFrom, pointTo);
                    if (this.numberOfLinesPerPair[pairStr] > 1 || linesPerPairCounter[pairStr] > 1) {
                        midPoint = this.findPerpendicularPointOnCurve(linesPerPairCounter[pairStr],
                            pointFrom.x, pointFrom.y, midPoint.x, midPoint.y)
                    }
                    line.quadraticCurveTo(midPoint.x, midPoint.y, bX, bY);
                    if (this.lines.drawArrows?.getBit(i))
                        this.canvasArrow(line, bX, bY, midPoint.x, midPoint.y);
                    this.ctx.beginPath();
                    this.ctx.stroke(line);
                    this.ctx.closePath();
                }
            }
        }
    }

    countLinesPerPair(): void {
        for (let i = 0; i < this.lines.from.length; i++) {
            const pairStr = this.getStringPair(this.lines.from[i], this.lines.to[i]);
            !this.numberOfLinesPerPair[pairStr] ? this.numberOfLinesPerPair[pairStr] = 1 : this.numberOfLinesPerPair[pairStr]++;
        }
    }

    getStringPair(point1: number, point2: number): string {
        let p1 = point1;
        let p2 = point2;
        if (point1 > point2) {
            p1 = point2;
            p2 = point1;
        }
        return `${p1}#${p2}`;
    }


    checkCoordsOnLine(x: number, y: number): number | null {
        for (let i = 0; i < this.paths.length; i++) {
            this.ctx.lineWidth = 5;
            const inStroke = this.ctx.isPointInStroke(this.paths[i], x * window.devicePixelRatio, y * window.devicePixelRatio);
            if (inStroke)
                return i;
        }
        return null;
    }

    midPoint(p1: DG.Point, p2: DG.Point): DG.Point {
        return new DG.Point(p1.x + (p2.x - p1.x) / 2, p1.y + (p2.y - p1.y) / 2);
    }

    findPerpendicularPointOnCurve(idx: number, x1: number, y1: number, x2: number, y2: number): DG.Point {
        let dx = x2 - x1;
        let dy = y2 - y1;
        const dist = Math.sqrt(dx * dx + dy * dy);
        dx /= dist;
        dy /= dist;
        const perpendicularLen = 50 * Math.ceil(idx / 2);
        return idx % 2 === 0 ?
            new DG.Point(x2 + (perpendicularLen / 2) * dy, y2 - (perpendicularLen / 2) * dx) :
            new DG.Point(x2 - (perpendicularLen / 2) * dy, y2 + (perpendicularLen / 2) * dx);
    }

    canvasArrow(path: Path2D, arrowEndX: number, arrowEndY: number, quadX: number, quadY: number): void {
        const arrowAngle = Math.atan2(quadX - arrowEndX, quadY - arrowEndY) + Math.PI;
        const arrowWidth = 15;
        path.moveTo(arrowEndX - (arrowWidth * Math.sin(arrowAngle - Math.PI / 10)),
            arrowEndY - (arrowWidth * Math.cos(arrowAngle - Math.PI / 10)));
        path.lineTo(arrowEndX, arrowEndY);
        path.lineTo(arrowEndX - (arrowWidth * Math.sin(arrowAngle + Math.PI / 10)),
            arrowEndY - (arrowWidth * Math.cos(arrowAngle + Math.PI / 10)));
    }

}