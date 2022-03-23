import { colorsForSurvivalChart } from "../constants/constants";
import { divideTo4Quadrants, divideTo6Quadrants, drawCensored, drawHorizontalInterval, numOfPatientsInQuadrant, renderLegend, renderNormalRanges, renderQuadrantsDefinitions, renderText, setCanvasStyle } from "./utils";
import * as DG from "datagrok-api/dg";


export function HysLawRenderLines(sp, xLimitAltAst, xLimitAp, yLimit, APColumnName) {
    let ctx = sp.getInfo()[ 'canvas' ].getContext('2d');
    setCanvasStyle(ctx, 'black', 0.2, '10px verdana');
    let xLimit = sp.getOptions().look.xColumnName === APColumnName ? xLimitAp : xLimitAltAst;
    let pointXLimit = sp.worldToScreen(xLimit, 0);
    let pointYLimit = sp.worldToScreen(0, yLimit);
    let left = sp.viewBox.x + 3;
    let right = sp.viewBox.x + sp.viewBox.width;
    let top = sp.viewBox.y;
    let bottom = sp.viewBox.y + sp.viewBox.height;
    renderText(ctx, `Normal(${numOfPatientsInQuadrant(sp, `${sp.getOptions().look.xColumnName} < ${xLimit} and ${sp.getOptions().look.yColumnName} < ${yLimit}`, sp.dataFrame.rowCount)})`,
        'bottom', 'left', left, bottom, left < pointXLimit.x && bottom > pointYLimit.y, 'blue')
    renderText(ctx, `Temple\'s Corollary(${numOfPatientsInQuadrant(sp, `${sp.getOptions().look.xColumnName} > ${xLimit} and ${sp.getOptions().look.yColumnName} < ${yLimit}`, sp.dataFrame.rowCount)})`,
        'bottom', 'right', right, bottom, right > pointXLimit.x && bottom > pointYLimit.y, 'blue')
    renderText(ctx, `Possible Hy\'s Law(${numOfPatientsInQuadrant(sp, `${sp.getOptions().look.xColumnName} > ${xLimit} and ${sp.getOptions().look.yColumnName} > ${yLimit}`, sp.dataFrame.rowCount)})`,
        'top', 'right', right, top, right > pointXLimit.x && top < pointYLimit.y, 'blue')
    renderText(ctx, `Hyperbilirubinemia(${numOfPatientsInQuadrant(sp, `${sp.getOptions().look.xColumnName} < ${xLimit} and ${sp.getOptions().look.yColumnName} > ${yLimit}`, sp.dataFrame.rowCount)})`,
        'top', 'left', left, top, left < pointXLimit.x && top < pointYLimit.y, 'blue')
    divideTo4Quadrants(ctx, sp, pointXLimit.x, pointYLimit.y);
}


export function BaselineEndpointRenderLines(sp, lowLevel, highLevel) {
    let ctx = sp.getInfo()[ 'canvas' ].getContext('2d');
    setCanvasStyle(ctx, 'black', 0.2, '10px verdana');
    let pointXHigh = sp.worldToScreen(highLevel, 0);
    let pointXLow = sp.worldToScreen(lowLevel, 0);
    let pointYHigh = sp.worldToScreen(0, highLevel);
    let pointYLow = sp.worldToScreen(0, lowLevel);
    let left = sp.viewBox.x + 3;
    let right = sp.viewBox.x + sp.viewBox.width;
    let top = sp.viewBox.y;
    let bottom = sp.viewBox.y + sp.viewBox.height;
    renderNormalRanges(ctx, left, right, top, bottom, pointXHigh, pointXLow, pointYHigh, pointYLow)
    renderQuadrantsDefinitions(ctx, left, right, top, bottom, pointXHigh, pointXLow, pointYHigh, pointYLow)
    divideTo6Quadrants(ctx, left, right, top, bottom, pointXHigh.x, pointXLow.x, pointYHigh.y, pointYLow.y);
}


export function KaplanMeierRenderLines(subjIDColumn, xColumn, yColumn, groupColumn, sp) {
    let ctx = sp.getInfo()[ 'canvas' ].getContext('2d');
    let groupsByColor = [];
    let groups = sp.dataFrame.groupBy([ groupColumn ])
        .aggregate();
    for (let i = 0; i < groups.rowCount; i++) {
        let group = groups.get(groupColumn, i)
        if (group) {
            let censored = [];
            let condition = groupColumn + ' = ' + group
            let groupRaws = sp.dataFrame.groupBy([ subjIDColumn, xColumn, yColumn, groupColumn ])
                .where(condition)
                .aggregate();
            ctx.beginPath();
            ctx.strokeStyle = colorsForSurvivalChart[ i ];
            let startPoint = sp.worldToScreen(0, 1);
            ctx.moveTo(startPoint.x, startPoint.y);
            for (let j = 0; j < groupRaws.rowCount; j++) {
                let survival = groupRaws.get(yColumn, j)
                let time = groupRaws.get(xColumn, j)
                if (j === 0) { groupsByColor.push({ group: group, color: colorsForSurvivalChart[ i ] }) }
                if (survival !== DG.FLOAT_NULL) {
                    let eventPoint = drawHorizontalInterval(sp, ctx, time, survival, startPoint.y);
                    ctx.lineTo(eventPoint.x, eventPoint.y);
                    startPoint = eventPoint
                } else {
                    censored.push({ x: time, y: startPoint.y })
                    if (j === groupRaws.rowCount - 1) {
                        drawHorizontalInterval(sp, ctx, time, survival, startPoint.y);
                    }
                }
                ctx.stroke();
            }
            drawCensored(ctx, sp, censored);
        }
    }
    renderLegend(sp, ctx, groupsByColor);
}


