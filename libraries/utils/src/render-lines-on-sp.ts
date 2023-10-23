import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import BitArray from './bit-array';

type MarkerSize = {
    sizeFrom: number,
    sizeTo: number
}

export type ILineSeries = {
    from: Uint32Array;
    to: Uint32Array;
    color?: string;             // common color. Use [colors] if you need individual colors per line
    colors?: string[];          // line colors. If they are the same for the series, use [color] instead.
    width?: number;             // common width. Use [colors] if you need individual widths per line
    widths?: Float32Array;      // line widths. If they are the same for the series, use [width] instead.
    opacity?: number;           // common opacity. Use [opacities] if you need individual opacities per line
    opacities?: Float32Array;   // line opacities. If they are the same for the series, use [opacity] instead
    drawArrows?: boolean;       // common parameter to draw arrows. Use [drawArrowsArr] if you need to draw arrows not for each line
    drawArrowsArr?: BitArray;   // individual parameter for each line. If they are the same for the series, use [drawArrows] instead

}

export class ScatterPlotLinesRenderer {
    sp: DG.ScatterPlotViewer;
    xAxisCol: DG.Column;
    yAxisCol: DG.Column;
    currentLineIdx = -1;
    lines: ILineSeries;
    lineClicked = new Subject<number>();
    lineHover = new Subject<number>();
    canvas: HTMLCanvasElement;
    ctx: CanvasRenderingContext2D;
    mouseOverLineId: number | null = null;
    multipleLinesCounts: Uint8Array;

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
        this.multipleLinesCounts = new Uint8Array(this.lines.from.length);

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
        const spLook = this.sp.getOptions().look;
        const individualLineStyles = this.lines.colors || this.lines.width || this.lines.opacities || this.lines.drawArrowsArr;
        if (!individualLineStyles) {
            this.ctx.lineWidth = this.lines.width ?? 1;
            this.ctx.strokeStyle = `rgba(${this.lines.color ?? '0,128,0'},${this.lines.opacity ?? 1})`;
        }
        const markerSizeCol = spLook['sizeColumnName'] ? this.sp.dataFrame.col(spLook['sizeColumnName']) : null;
        for (let i = 0; i < this.lines.from.length; i++) {
            const { sizeFrom, sizeTo } = this.getMarkersSizes(spLook, markerSizeCol, i);
            const pointFrom = this.sp.worldToScreen(this.xAxisCol.get(this.lines.from[i]), this.yAxisCol.get(this.lines.from[i]));
            let aX = pointFrom?.x;
            let aY = pointFrom?.y;
            const pointTo = this.sp.worldToScreen(this.xAxisCol.get(this.lines.to[i]), this.yAxisCol.get(this.lines.to[i]));
            let bX = pointTo?.x;
            let bY = pointTo?.y;
            this.ctx.beginPath();
            if (aX && aY && bX && bY) {
                if (individualLineStyles) {
                    const color = this.lines.colors?.[i] ? this.lines.colors?.[i] : '0,128,0';
                    const opacity = this.lines.opacities?.[i] ? this.lines.opacities?.[i] : 1;
                    this.ctx.strokeStyle = `rgba(${color},${opacity})`;
                    this.ctx.lineWidth = this.lines.widths?.[i] ? this.lines.widths?.[i] : 1;
                }
                const multiLines = this.multipleLinesCounts[i];
                let controlPoint: DG.Point | null = null;
                if (multiLines) {
                    const startPointWithMarker = this.getPointOnDistance(aX, aY, bX, bY, sizeTo);
                    const endtPointWithMarker = this.getPointOnDistance(bX, bY, aX, aY, sizeFrom);
                    aX = startPointWithMarker.x;
                    aY = startPointWithMarker.y;
                    bX = endtPointWithMarker.x;
                    bY = endtPointWithMarker.y;
                    controlPoint = this.lines.from[i] > this.lines.to[i] ?
                        this.findControlPoint(multiLines, aX, aY, bX, bY, i) :
                        this.findControlPoint(multiLines, bX, bY, aX, aY, i);
                    this.ctx.moveTo(aX!, aY!);
                    this.ctx.quadraticCurveTo(controlPoint.x, controlPoint.y, bX, bY);
                } else {
                    this.ctx.moveTo(aX!, aY!);
                    this.ctx.lineTo(bX, bY);
                }
                if (this.lines.drawArrows ?? this.lines.drawArrowsArr?.getBit(i)) {
                    const arrowPoint = !multiLines ? this.getPointOnDistance(aX, aY, bX, bY, sizeTo) : null;
                    const arrowCPX = multiLines ? controlPoint!.x : aX;
                    const arrowCPY = multiLines ? controlPoint!.y : aY;
                    this.canvasArrow(this.ctx, arrowPoint?.x ?? bX, arrowPoint?.y ?? bY, arrowCPX, arrowCPY);
                }
                this.ctx.stroke();
                this.ctx.closePath();
            }
        }
        this.fillLeftBottomRect();
    }


    getMarkersSizes(spLook: any, markerSizeCol: DG.Column | null, i: number): MarkerSize {
        let sizeFrom = 3;
        let sizeTo = 3;
        if (markerSizeCol) {
            sizeFrom = (spLook.markerMinSize + (spLook.markerMaxSize - spLook.markerMinSize) * markerSizeCol.scale(this.lines.from[i])) / 2;
            sizeTo = (spLook.markerMinSize + (spLook.markerMaxSize - spLook.markerMinSize) * markerSizeCol.scale(this.lines.to[i])) / 2;
        } else if (spLook.markerDefaultSize) {
            sizeFrom = spLook.markerDefaultSize / 2;
            sizeTo = spLook.markerDefaultSize / 2;
        }
        return { sizeFrom, sizeTo };
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
        const arrayIdxsBitArray = new BitArray(this.lines.from.length);
        arrayIdxsBitArray.setAll(true);
        for (let i = -1; (i = arrayIdxsBitArray.findNext(i)) !== -1;) {
            let firstLineIdx = i;
            let p1 = this.lines.from[firstLineIdx];
            let p2 = this.lines.to[firstLineIdx];
            let linesPerPair = 1;
            for (let j = i; (j = arrayIdxsBitArray.findNext(j)) !== -1;) {
                const pointToCompare1 = this.lines.from[j];
                const pointToCompare2 = this.lines.to[j];
                if (pointToCompare1 === p1 && pointToCompare2 === p2 ||
                    pointToCompare2 === p1 && pointToCompare1 === p2) {
                    this.multipleLinesCounts[j] = ++linesPerPair;
                    arrayIdxsBitArray.setBit(j, false, false);
                }
            }
            if (linesPerPair > 1)
                this.multipleLinesCounts[firstLineIdx] = 1;
            arrayIdxsBitArray.setBit(i, false, false);
        }
    }

    checkCoordsOnLine(x: number, y: number): number | null {
        let candidateIdx = null;
        let minDist = null;
        let dist = null;
        const spLook = this.sp.getOptions().look;
        const markerSizeCol = spLook['sizeColumnName'] ? this.sp.dataFrame.col(spLook['sizeColumnName']) : null;
        for (let i = 0; i < this.lines.from.length; i++) {
            const { sizeFrom, sizeTo } = this.getMarkersSizes(spLook, markerSizeCol, i);
            const pFrom = this.sp.worldToScreen(this.xAxisCol.get(this.lines.from[i]), this.yAxisCol.get(this.lines.from[i]));
            const pTo = this.sp.worldToScreen(this.xAxisCol.get(this.lines.to[i]), this.yAxisCol.get(this.lines.to[i]));
            if (this.multipleLinesCounts[i]) {
                const fromMarker = this.getPointOnDistance(pFrom.x, pFrom.y, pTo.x, pTo.y, sizeTo);
                const toMarker = this.getPointOnDistance(pTo.x, pTo.y, pFrom?.x, pFrom?.y, sizeFrom);
                const controlPoint = this.lines.from[i] > this.lines.to[i] ?
                    this.findControlPoint(this.multipleLinesCounts[i], fromMarker.x, fromMarker.y, toMarker.x, toMarker.y, i) :
                    this.findControlPoint(this.multipleLinesCounts[i], toMarker.x, toMarker.y, fromMarker.x, fromMarker.y, i);
                dist = this.calculateDistToCurveLine(i, x, y, fromMarker, toMarker, controlPoint);
            } else {
                dist = this.calculateDistToStraightLine(x, y, pFrom, pTo);
            }
            if ((!minDist && dist !== null && dist < 5) || minDist && dist !== null && dist < minDist) {
                minDist = dist;
                candidateIdx = i;
            }
        }
        return candidateIdx;
    }

    calculateDistToStraightLine(x: number, y: number, p1: DG.Point, p2: DG.Point): number | null {
        /* calculating coordinates of a rect around a line. If cursor coords are outside this rect - assume that
        point is not on line and do not calculate distance to line */
        let xMin = Math.min(p1.x, p2.x);
        let xMax = Math.max(p1.x, p2.x);
        let yMin = Math.min(p1.y, p2.y);
        let yMax = Math.max(p1.y, p2.y);

        //adding a couple of pixels to increase the width/height of the rect
        const threshold = 2;
        return x >= xMin - threshold && x <= xMax + threshold && y >= yMin - threshold && y <= yMax + threshold ?
            Math.abs(Math.hypot(p1.x - x, p1.y - y) + Math.hypot(p2.x - x, p2.y - y) - Math.hypot(p1.x - p2.x, p1.y - p2.y)) :
            null;
    }

    calculateDistToCurveLine(i: number, x: number, y: number, p1: DG.Point, p2: DG.Point,
        pc: DG.Point): number | null {
        /* calculating coordinates of a rect around a line. If cursor coords are outside this shape - assume that
        point is not on line and do not calculate distance to line */
        let xMin = Math.min(p1.x, p2.x, pc.x);
        let xMax = Math.max(p1.x, p2.x, pc.x);
        let yMin = Math.min(p1.y, p2.y, pc.y);
        let yMax = Math.max(p1.y, p2.y, pc.y);

        //adding a couple of pixels to increase the width/height of the rect
        const threshold = 2;
        if (x >= xMin - threshold && x <= xMax + threshold && y >= yMin - threshold && y <= yMax + threshold) {
            const w = xMax - xMin;
            const h = yMax - yMin;
            return this.calculateDistToCurveInRect(x, y, p1, pc, p2, w, h);
        }
        return null;
    }

    calculateDistToCurveInRect(x: number, y: number, p0: DG.Point, p1: DG.Point, p2: DG.Point,
        w: number, h: number): number {
        const stepLen = 3;
        const stepsNum = Math.floor((w + h) / stepLen);
        const deltaT = 1 / stepsNum;
        const arrX = new Uint32Array(stepsNum);
        const arrY = new Uint32Array(stepsNum);
        let maxHW = new Uint32Array(stepsNum);
        let minSumHW: number | null = null;
        const candidateIdxs = new BitArray(stepsNum);
        for (let i = 0; i < arrX.length; i++) {
            const t = i * deltaT;
            const xOnCurve = Math.pow((1 - t), 2) * p0.x + 2 * t * (1 - t) * p1.x + Math.pow(t, 2) * p2.x;
            const yOnCurve = Math.pow((1 - t), 2) * p0.y + 2 * t * (1 - t) * p1.y + Math.pow(t, 2) * p2.y;
            const rectW = Math.abs(x - xOnCurve);
            const rectH = Math.abs(y - yOnCurve);
            const sumHW = rectW + rectH;
            if (!minSumHW || minSumHW > sumHW)
                minSumHW = sumHW;
            maxHW[i] = Math.max(rectW, rectH);
            arrX[i] = xOnCurve;
            arrY[i] = yOnCurve;
        }
        for (let i = 0; i < arrX.length; i++) {
            if (maxHW[i] < minSumHW!)
                candidateIdxs.setBit(i, true, false);
        }
        let minDist: number | null = null;
        for (let j = -1; (j = candidateIdxs.findNext(j)) !== -1;) {
            const dist = Math.hypot((arrX[j] - x), (arrY[j] - y));
            if (!minDist || minDist > dist)
                minDist = dist;
        }
        return minDist!;
    }

    getPointOnDistance(p1x: number, p1y: number, p2x: number, p2y: number, distance: number): DG.Point {
        const p1p2d = Math.sqrt(Math.pow(p2x - p1x, 2) + Math.pow(p2y - p1y, 2));
        const dx = (p2x - p1x) / p1p2d;
        const dy = (p2y - p1y) / p1p2d;
        const p3x = p2x - distance * dx;
        const p3y = p2y - distance * dy;

        return new DG.Point(p3x, p3y);
    }

    findControlPoint(idx: number, x1: number, y1: number, x2: number, y2: number, i?: number): DG.Point {
        const midX = x1 + (x2 - x1) / 2;
        const midY = y1 + (y2 - y1) / 2;
        let dx = midX - x1;
        let dy = midY - y1;
        const dist = Math.sqrt(dx * dx + dy * dy);
        dx /= dist;
        dy /= dist;
        const perpendicularLen = 50 * Math.ceil(idx / 2);
        return idx % 2 === 0 ?
            new DG.Point(midX + (perpendicularLen / 2) * dy, midY - (perpendicularLen / 2) * dx) :
            new DG.Point(midX - (perpendicularLen / 2) * dy, midY + (perpendicularLen / 2) * dx);
    }

    canvasArrow(path: CanvasRenderingContext2D, arrowEndX: number, arrowEndY: number, quadX: number, quadY: number): void {
        const arrowAngle = Math.atan2(quadX - arrowEndX, quadY - arrowEndY) + Math.PI;
        const arrowWidth = 15;
        path.moveTo(arrowEndX - (arrowWidth * Math.sin(arrowAngle - Math.PI / 10)),
            arrowEndY - (arrowWidth * Math.cos(arrowAngle - Math.PI / 10)));
        path.lineTo(arrowEndX, arrowEndY);
        path.lineTo(arrowEndX - (arrowWidth * Math.sin(arrowAngle + Math.PI / 10)),
            arrowEndY - (arrowWidth * Math.cos(arrowAngle + Math.PI / 10)));
    }

}