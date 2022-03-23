import { SUBJECT_ID } from "../constants/columns-constants";

export function renderText(ctx, text, baseline, align, x, y, condition, color = 'black') {
  if (condition) {
    ctx.textBaseline = baseline;
    ctx.textAlign = align;
    ctx.fillStyle = color
    ctx.fillText(text, x, y);
  }
}

export function setCanvasStyle(ctx, color, lineWidth, font) {
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth;
  ctx.font = font;
}


export function divideTo4Quadrants(ctx, sp, limitX, limitY) {
  ctx.beginPath();
  ctx.moveTo(sp.viewBox.x, limitY);
  ctx.lineTo(sp.viewBox.x + sp.viewBox.width, limitY);
  ctx.moveTo(limitX, sp.viewBox.y);
  ctx.lineTo(limitX, sp.viewBox.y + sp.viewBox.height);
  ctx.stroke();
}

export function divideTo6Quadrants(ctx, left, right, top, bottom, xHigh, xLow, yHigh, yLow) {
  ctx.beginPath();
  ctx.moveTo(xLow, bottom);
  ctx.lineTo(xLow, top);
  ctx.moveTo(xHigh, bottom);
  ctx.lineTo(xHigh, top);
  ctx.moveTo(right, yLow);
  ctx.lineTo(left, yLow);
  ctx.moveTo(right, yHigh);
  ctx.lineTo(left, yHigh);
  ctx.stroke();
}


export function numOfPatientsInQuadrant(sp, condition, totalNum) {
  return sp.dataFrame.groupBy([ SUBJECT_ID ])
    .where(condition)
    .aggregate().rowCount;
}


export function renderNormalRanges(ctx, left, right, top, bottom, pointXHigh, pointXLow, pointYHigh, pointYLow) {
  renderText(ctx, 'LLNR', 'bottom', 'left', left, pointYLow.y, true);
  renderText(ctx, 'ULNR', 'bottom', 'left', left, pointYHigh.y, true);
  renderText(ctx, 'LLNR', 'top', 'left', pointXLow.x, top, true);
  renderText(ctx, 'ULNR', 'top', 'left', pointXHigh.x, top, true);
}

export function renderQuadrantsDefinitions(ctx, left, right, top, bottom, pointXHigh, pointXLow, pointYHigh, pointYLow) {
  renderText(ctx, 'Low-Low', 'top', 'right', pointXLow.x, pointYLow.y, true, 'blue');
  renderText(ctx, 'Low-Normal', 'top', 'right', pointXLow.x, pointYHigh.y, true, 'blue');
  renderText(ctx, 'Low-High', 'top', 'right', pointXLow.x, pointYHigh.y - top > 10 ? top : pointYHigh.y - 10, true, 'blue');
  renderText(ctx, 'Normal-Low', 'top', 'right', pointXHigh.x, pointYLow.y, true, 'blue');
  renderText(ctx, 'Normal-Normal', 'top', 'right', pointXHigh.x, pointYHigh.y, true, 'blue');
  renderText(ctx, 'Normaml-High', 'top', 'right', pointXHigh.x, pointYHigh.y - top > 10 ? top : pointYHigh.y - 10, true, 'blue');
  renderText(ctx, 'High-Low', 'top', 'right', right - pointXHigh.x > 60 ? right : pointXHigh, pointYLow.y, true, 'blue');
  renderText(ctx, 'High-Normal', 'top', 'right', right - pointXHigh.x > 60 ? right : pointXHigh, pointYHigh.y, true, 'blue');
  renderText(ctx, 'High-High', 'top', 'right', right - pointXHigh.x > 60 ? right : pointXHigh, pointYHigh.y - top > 10 ? top : pointYHigh.y - 10, true, 'blue');
}


export function drawHorizontalInterval(sp, ctx, time, survival, startY) {
  const eventPoint = sp.worldToScreen(time, survival);
  ctx.lineTo(eventPoint.x, startY);
  return eventPoint;
}

export function drawCensored(ctx, sp, censored) {
  censored.forEach(item => {
    ctx.beginPath();
    const censoredPoint = sp.worldToScreen(item.x, 0);
    ctx.moveTo(censoredPoint.x - 4, item.y + 2);
    ctx.lineTo(censoredPoint.x + 4, item.y + 2);
    ctx.moveTo(censoredPoint.x, item.y - 2);
    ctx.lineTo(censoredPoint.x, item.y + 6);
    ctx.stroke();
  })
}


export function renderLegend(sp, ctx, groupsByColor) {
  let left = sp.viewBox.x + 3;
  let bottom = sp.viewBox.y + sp.viewBox.height;
  renderText(ctx, '+ \u2013 censored', 'bottom', 'left', left, bottom, 'black');
  groupsByColor.forEach((item, index) => {
    renderText(ctx, item.group, 'bottom', 'left', left, bottom - (index + 1) * 20, true, item.color)
  })
}
