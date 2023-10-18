import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import BitArray from './bit-array';

export type ILineSeries = {
    from: Uint32Array;
    to: Uint32Array;
    color?: string;             // common color. Use [colors] if you need individual colors per line
    colors?: string[];          // line colors. If they are the same for the series, use [color] instead.
    width?: number;             // common width. Use [colors] if you need individual widths per line
    widths?: Float32Array;      // line widths. If they are the same for the series, use [width] instead.
    opacity?: number;           // common opacity. Use [opacities] if you need individual opacities per line
    opacities?: Float32Array;   // line opacities. If they are the same for the series, use [opacity] instead
    drawArrows?: Boolean;       // common parameter to draw arrows. Use [drawArrowsArr] if you need to draw arrows not for each line
    drawArrowsArr?: BitArray;   // individual parameter for each line. If they are the same for the series, use [drawArrows] instead

}

export class ScatterPlotLinesRenderer {
    sp: DG.ScatterPlotViewer;
    xAxisCol: DG.Column;
    yAxisCol: DG.Column;;
    currentLineIdx = -1;
    lines: ILineSeries;
    lineClicked = new Subject<number>();
    lineHover = new Subject<number>();
    canvas: HTMLCanvasElement;
    ctx: CanvasRenderingContext2D;
    paths: Path2D[];
    mouseOverLineId: number | null = null;
    multipleLinesIndices: Uint32Array;
    changePointsOrder: BitArray;

    get currentLineId(): number {
        return this.currentLineIdx;
    }

    set currentLineId(id: number) {
        this.currentLineIdx = id;
        this.sp.render(this.canvas.getContext('2d') as CanvasRenderingContext2D);
    }

    constructor(sp: DG.ScatterPlotViewer, xAxis: string, yAxis: string, lines: ILineSeries) {
        this.sp = sp;
        this.xAxisCol = this.sp.dataFrame!.columns.byName(xAxis);
        this.yAxisCol = this.sp.dataFrame!.columns.byName(yAxis);
        this.lines = lines;
        this.canvas = this.sp.getInfo()['canvas'];
        this.ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        this.paths = new Array<Path2D>(this.lines.from.length);
        this.multipleLinesIndices = new Uint32Array(this.lines.from.length);
        this.changePointsOrder = new BitArray(this.lines.from.length, false);

        this.createMultiLinesIndices();

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
        for (let i = 0; i < this.lines.from.length; i++) {
            const pointsOrderChanged = this.changePointsOrder.getBit(i);
            let pointFromIdx = pointsOrderChanged ? this.lines.to[i] : this.lines.from[i];
            let pointToIdx = pointsOrderChanged ? this.lines.from[i] : this.lines.to[i];
            const pointFrom = this.sp.worldToScreen(this.xAxisCol.get(pointFromIdx), this.yAxisCol.get(pointFromIdx));
            const aX = pointFrom?.x;
            const aY = pointFrom?.y;
            const pointTo = this.sp.worldToScreen(this.xAxisCol.get(pointToIdx), this.yAxisCol.get(pointToIdx));
            const bX = pointTo?.x;
            const bY = pointTo?.y;
            if (aX && aY && bX && bY) {
                const line = new Path2D();
                this.paths[i] = line;
                line.moveTo(aX!, aY!);
                const color = this.lines.color ?? this.lines.colors?.[i] ? this.lines.colors?.[i] : '0,128,0';
                const opacity = this.lines.opacity ?? this.lines.opacities?.[i] ? this.lines.opacities?.[i] : 1;
                this.ctx.strokeStyle = `rgba(${color},${opacity})`;
                this.ctx.lineWidth = this.lines.width ? this.lines.width : this.lines.widths?.[i] ? i === this.currentLineIdx ? this.lines.widths?.[i] + 2 :
                    this.lines.widths?.[i] : i === this.currentLineIdx ? 3 : 1;
                let midPoint: DG.Point = this.midPoint(pointFrom, pointTo);
                if (this.multipleLinesIndices[i]) {
                    midPoint = this.findPerpendicularPointOnCurve(this.multipleLinesIndices[i],
                        pointFrom.x, pointFrom.y, midPoint.x, midPoint.y)
                }
                line.quadraticCurveTo(midPoint.x, midPoint.y, bX, bY);
                if (this.lines.drawArrows ?? this.lines.drawArrowsArr?.getBit(i))
                    this.canvasArrow(line, pointsOrderChanged ? aX : bX, pointsOrderChanged ? aY : bY, midPoint.x, midPoint.y);
                this.ctx.beginPath();
                this.ctx.stroke(line);
                this.ctx.closePath();
            }
        }
        this.fillLeftBottomRect();
    }

    fillLeftBottomRect() {
        const rect = new Path2D();
        rect.rect(this.sp.yAxisBox.minX, this.sp.yAxisBox.maxY, this.sp.yAxisBox.width, this.sp.xAxisBox.height);
        this.ctx.fillStyle = `white`;
        this.ctx.beginPath();
        this.ctx.fill(rect);
        this.ctx.closePath();
    }

    createMultiLinesIndices(): void {
        const arrayIdxs = new Array<number>(this.lines.from.length);
        for (let i = 0; i < arrayIdxs.length; i++)
            arrayIdxs[i] = i;
        while (arrayIdxs.length) {
            let firstLineIdx = arrayIdxs[0];
            let p1 = this.lines.from[firstLineIdx];
            let p2 = this.lines.to[firstLineIdx];
            if (p1 > p2)
                this.changePointsOrder.setBit(firstLineIdx, true, false);
            let linesPerPair = 1;
            for (let i = 1; i < arrayIdxs.length; i++) {
                const pointToCompare1 = this.lines.from[arrayIdxs[i]];
                const pointToCompare2 = this.lines.to[arrayIdxs[i]];
                if (pointToCompare1 === p1 && pointToCompare2 === p2 ||
                    pointToCompare2 === p1 && pointToCompare1 === p2) {
                    if (pointToCompare1 > pointToCompare2)
                        this.changePointsOrder.setBit(arrayIdxs[i], true);
                    this.multipleLinesIndices[arrayIdxs[i]] = ++linesPerPair;
                    arrayIdxs.splice(i, 1);
                    i--;
                }
            }
            if (linesPerPair > 1)
                this.multipleLinesIndices[firstLineIdx] = 1;
            arrayIdxs.splice(0, 1);
        }
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