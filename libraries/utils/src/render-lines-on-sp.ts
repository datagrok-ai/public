import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';

export enum LineDirection {
    lower = 'lower',
    higher = 'higher',
}

export type ILineOpts = {
    id: number;
    color?: string;
    directionCol?: string
    direction?: LineDirection;
    opacity?: number;
    width?: number;
}

export type IConnectedPoints = { [key: number]: ILineOpts[] }

export type ILines = { [key: number]: IConnectedPoints }

export class ScatterPlotLinesRenderer {
    sp: DG.ScatterPlotViewer;
    xAxis: string;
    yAxis: string;
    currentLineIdx = -1;
    lines: ILines;
    lineClicked = new Subject<number>();
    lineHover = new Subject<number>();
    canvas: any;
    paths: { [key: number]: Path2D } = {};
    mouseOverLineId: number | null = null;

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
        this.canvas = (this.sp.getInfo() as { [index: string]: any })['canvas'];

        this.canvas.addEventListener('mousedown', (event: MouseEvent) => {
            if (this.mouseOverLineId !== null)
                this.lineClicked.next(this.mouseOverLineId);
        });

        this.canvas.addEventListener('mousemove', (event: MouseEvent) => {
            this.mouseOverLineId = this.checkCursorOnLine(event.offsetX, event.offsetY);
            if (this.mouseOverLineId !== null)
                this.lineHover.next(this.mouseOverLineId);
        });

        sp.onEvent('d4-before-draw-scene')
            .subscribe((_: any) => {
                this.renderLines();
            });

    }

    renderLines(): void {
        const ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        const x = this.sp.dataFrame!.columns.byName(this.xAxis);
        const y = this.sp.dataFrame!.columns.byName(this.yAxis);
        for (let p1 of Object.keys(this.lines).map(Number)) {

            if (this.sp.dataFrame.filter.get(p1)) {
                const pointFrom = this.sp.worldToScreen(x.get(p1), y.get(p1));
                const aX = pointFrom?.x;
                const aY = pointFrom?.y;
                for (let p2 of Object.keys(this.lines[p1]).map(Number)) {
                    if (this.sp.dataFrame.filter.get(p2)) {
                        const pointTo = this.sp.worldToScreen(x.get(p2), y.get(p2));
                        const bX = pointTo?.x;
                        const bY = pointTo?.y;
                        if (aX && aY && bX && bY) {
                            for (let lineIdx = 0; lineIdx < this.lines[p1][p2].length; lineIdx++) {
                                const lineOpts = this.lines[p1][p2][lineIdx];
                                const line = new Path2D();
                                this.paths[lineOpts.id] = line;
                                line.moveTo(aX!, aY!);
                                const color = lineOpts.color ?? '0,128,0';
                                const opacity = lineOpts.opacity ?? 1;
                                ctx.strokeStyle = `rgba(${color},${opacity})`;
                                ctx.lineWidth = lineOpts.width ? lineOpts.id === this.currentLineIdx ? lineOpts.width + 2 : lineOpts.width :
                                    lineOpts.id === this.currentLineIdx ? 3 : 1;
                                let midPoint: DG.Point = this.midPointBtw(pointFrom, pointTo);
                                if (lineIdx > 0 || this.lines[p1][p2].length > 1) {
                                    midPoint = this.findPerpendicularPointOnCurve(lineIdx,
                                        pointFrom.x, pointFrom.y, midPoint.x, midPoint.y)
                                }
                                line.quadraticCurveTo(midPoint.x, midPoint.y, bX, bY);
                                if (lineOpts.directionCol) {
                                    const arrowPoint = this.getArrowPoint(lineOpts, p1, p2, pointFrom, pointTo);
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
            }
        }
    }

    checkCursorOnLine(x: number, y: number): number | null {
        const ctx = this.canvas.getContext('2d') as CanvasRenderingContext2D;
        for (let id of Object.keys(this.paths).map(Number)) {
            ctx.lineWidth = 5;
            const inStroke = ctx.isPointInStroke(this.paths[id], x * window.devicePixelRatio, y * window.devicePixelRatio);
            if (inStroke)
                return id;
        }
        return null;
    }

    midPointBtw(p1: DG.Point, p2: DG.Point): DG.Point {
        return new DG.Point(p1.x + (p2.x - p1.x) / 2, p1.y + (p2.y - p1.y) / 2);
    }

    findPerpendicularPointOnCurve(idx: number, x1: number, y1: number, x2: number, y2: number) {
        let dx = x2 - x1
        let dy = y2 - y1
        const dist = Math.sqrt(dx * dx + dy * dy)
        dx /= dist
        dy /= dist
        const perpendicularLen = 50 * Math.floor(idx / 2) + 50;
        return idx % 2 === 0 ?
            new DG.Point(x2 + (perpendicularLen / 2) * dy, y2 - (perpendicularLen / 2) * dx) :
            new DG.Point(x2 - (perpendicularLen / 2) * dy, y2 + (perpendicularLen / 2) * dx);
    }

    canvasArrow(path: Path2D, arrowEndX: number, arrowEndY: number, quadX: number, quadY: number) {
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
            const direction = lineOpts.direction ?? LineDirection.higher;
            switch (direction) {
                case LineDirection.higher:
                    return compareCol.get(p1) > compareCol.get(p2) ? pointFrom : pointTo;
                case LineDirection.lower:
                    return compareCol.get(p1) > compareCol.get(p2) ? pointTo : pointFrom;
                default:
                    return null;
            }
        }
        return null;
    }

}