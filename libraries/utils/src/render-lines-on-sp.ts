import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';

export enum LineDirection {
    pointToLowerValue = 'lower',
    pointToHigherValue = 'higher',
}

export type ILineOpts = {
    id: number;
    color?: string;
    directionCol?: string
    direction?: LineDirection;
    opacity?: number;
    width?: number;
}

export type ILines = {
    pairs: Uint32Array[];
    linesStyles: ILineOpts[];
}


export class ScatterPlotLinesRenderer {
    sp: DG.ScatterPlotViewer;
    xAxis: string;
    yAxis: string;
    currentLineIdx = -1;
    lines: ILines;
    lineClicked = new Subject<number>();
    lineHover = new Subject<number>();
    canvas: HTMLCanvasElement;
    paths: { [key: number]: Path2D } = {};
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
        this.xAxis = xAxis;
        this.yAxis = yAxis;
        this.lines = lines;
        this.canvas = this.sp.getInfo()['canvas'];
        this.countLinesPerPairAndSort();

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
        const ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        const x = this.sp.dataFrame!.columns.byName(this.xAxis);
        const y = this.sp.dataFrame!.columns.byName(this.yAxis);
        const linesPerPairCounter: { [key: string]: number } = {};
        for (let i = 0; i < this.lines.pairs.length; i++) {
            const pair = this.lines.pairs[i];
            const lineOpts = this.lines.linesStyles[i];
            if (this.sp.dataFrame.filter.get(pair[0]) && this.sp.dataFrame.filter.get(pair[1])) {
                const pointFrom = this.sp.worldToScreen(x.get(pair[0]), y.get(pair[0]));
                const aX = pointFrom?.x;
                const aY = pointFrom?.y;
                const pointTo = this.sp.worldToScreen(x.get(pair[1]), y.get(pair[1]));
                const bX = pointTo?.x;
                const bY = pointTo?.y;
                if (aX && aY && bX && bY) {
                    const pairStr = this.getStringPair(pair[0], pair[1]);
                    !linesPerPairCounter[pairStr] ? linesPerPairCounter[pairStr] = 1 : linesPerPairCounter[pairStr]++;
                    const line = new Path2D();
                    this.paths[lineOpts.id] = line;
                    line.moveTo(aX!, aY!);
                    const color = lineOpts.color ?? '0,128,0';
                    const opacity = lineOpts.opacity ?? 1;
                    ctx.strokeStyle = `rgba(${color},${opacity})`;
                    ctx.lineWidth = lineOpts.width ? lineOpts.id === this.currentLineIdx ? lineOpts.width + 2 : lineOpts.width :
                        lineOpts.id === this.currentLineIdx ? 3 : 1;
                    let midPoint: DG.Point = this.midPoint(pointFrom, pointTo);
                    if (this.numberOfLinesPerPair[pairStr] > 1 || linesPerPairCounter[pairStr] > 1) {
                        midPoint = this.findPerpendicularPointOnCurve(linesPerPairCounter[pairStr],
                            pointFrom.x, pointFrom.y, midPoint.x, midPoint.y)
                    }
                    line.quadraticCurveTo(midPoint.x, midPoint.y, bX, bY);
                    if (lineOpts.directionCol) {
                        const arrowPoint = this.getArrowPoint(lineOpts, pair[0], pair[1], pointFrom, pointTo);
                        if (arrowPoint)
                            this.canvasArrow(line, arrowPoint.x, arrowPoint.y, midPoint.x, midPoint.y);
                    }
                    ctx.beginPath();
                    ctx.stroke(line);
                    ctx.closePath();
                }
            }
        }

    }

    countLinesPerPairAndSort(): void {
        for (let i = 0; i < this.lines.pairs.length; i++) {
            const pair = this.lines.pairs[i];
            const p1 = pair[0] > pair[1] ? pair[1] : pair[0];
            const p2 = pair[0] > pair[1] ? pair[0] : pair[1];
            pair[0] = p1;
            pair[1] = p2;
            const pairStr = `${this.lines.pairs[i][0]}#${this.lines.pairs[i][1]}`;
            if (!this.numberOfLinesPerPair[pairStr])
                this.numberOfLinesPerPair[pairStr] = 1;
            else
                this.numberOfLinesPerPair[pairStr]++;
        }
    }

    getStringPair(point1: number, point2: number): string {
        return `${point1}#${point2}`;
    }


    checkCoordsOnLine(x: number, y: number): number | null {
        const ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        for (let id of Object.keys(this.paths).map(Number)) {
            ctx.lineWidth = 5;
            const inStroke = ctx.isPointInStroke(this.paths[id], x * window.devicePixelRatio, y * window.devicePixelRatio);
            if (inStroke)
                return id;
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
        const perpendicularLen = 50 * Math.floor(idx / 2) + 50;
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

    getArrowPoint(lineOpts: ILineOpts, p1: number, p2: number, pointFrom: DG.Point, pointTo: DG.Point): DG.Point | null {
        const compareCol = this.sp.dataFrame!.col(lineOpts.directionCol!);
        if (compareCol) {
            const direction = lineOpts.direction ?? LineDirection.pointToHigherValue;
            switch (direction) {
                case LineDirection.pointToHigherValue:
                    return compareCol.get(p1) > compareCol.get(p2) ? pointFrom : pointTo;
                case LineDirection.pointToLowerValue:
                    return compareCol.get(p1) > compareCol.get(p2) ? pointTo : pointFrom;
                default:
                    return null;
            }
        }
        return null;
    }

}