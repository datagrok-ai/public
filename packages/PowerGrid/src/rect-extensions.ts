import {Rect} from "datagrok-api/src/grid";

declare module "datagrok-api/src/grid" {
  interface Rect {
    getTop(height: number): Rect;

    getBottom(height: number): Rect;

    getLeft(width: number): Rect;

    getRight(width: number): Rect;

    getTopLeft(width: number, height: number): Rect;

    getTopRight(width: number, height: number): Rect;

    getBottomLeft(width: number, height: number): Rect;

    getBottomRight(width: number, height: number): Rect;

    //--

    cutLeft(dw: number): Rect;

    cutTop(dh: number): Rect;

    cutBottom(dh: number): Rect;

    cutRight(dw: number): Rect;


    // --

    below(height: number): Rect;

    above(height: number): Rect;

    // --

    toTheLeft(width: number): Rect;

    toTheRight(width: number): Rect;

    toTheTop(height: number): Rect;

    toTheBottom(height: number): Rect;

    // --

    getTopScaled(ratio: number): Rect;

    getBottomScaled(ratio: number): Rect;

    getLeftScaled(ratio: number): Rect;

    getRightScaled(ratio: number): Rect;

    // --

    getTopPart(count: number, index: number): Rect;

    getLeftPart(count: number, index: number): Rect;

    getGridPart(xCount: number, yCount: number, x: number, y: number): Rect;

    // --

    inflate(dx: number, dy: number): Rect;

    inflateSize(dw: number, dh: number): Rect;

    inflateRel(dxRatio: number, dyRatio: number): Rect;

    // --

    fitSquare(): Rect;
  }
}

Rect.prototype.getTop = function (height) {
  return new Rect(this.left, this.top, this.width, height);
}

Rect.prototype.getBottom = function (height) {
  return new Rect(this.left, this.bottom - height, this.width, height)
}

Rect.prototype.getLeft = function (width) {
  return new Rect(this.left, this.top, width, this.height);
}

Rect.prototype.getRight = function (width) {
  return new Rect(this.right - width, this.top, width, this.height);
}

Rect.prototype.getTopLeft = function (width, height) {
  return new Rect(this.left, this.top, width, height);
}

Rect.prototype.getTopRight = function (width, height) {
  return new Rect(this.right - width, this.top, width, height);
}

Rect.prototype.getBottomLeft = function (width, height) {
  return new Rect(this.left, this.bottom - height, width, height);
}

Rect.prototype.getBottomRight = function (width, height) {
  return new Rect(this.right - width, this.bottom - height, width, height);
}

Rect.prototype.cutLeft = function (dw) {
  return new Rect(this.left + dw, this.top, this.width - dw, this.height);
}

Rect.prototype.cutTop = function (dh) {
  return new Rect(this.left, this.top + dh, this.width, this.height - dh);
}

Rect.prototype.cutBottom = function (dh) {
  return new Rect(this.left, this.top, this.width, this.height - dh);
}

Rect.prototype.cutRight = function (dw) {
  return new Rect(this.left, this.top, this.width - dw, this.height);
}

// --

Rect.prototype.below = function (height) {
  return new Rect(this.left, this.bottom, this.width, height);
}

Rect.prototype.above = function (height) {
  return new Rect(this.left, this.top - height, this.width, height);
}

// --

Rect.prototype.toTheLeft = function (width) {
  return new Rect(this.left - width, this.top, width, this.height);
}

Rect.prototype.toTheRight = function (width) {
  return new Rect(this.right, this.top, width, this.height);
}

Rect.prototype.toTheTop = function (height) {
  return this.above(height);
}

Rect.prototype.toTheBottom = function (height) {
  return this.below(height);
}

//--

Rect.prototype.getTopScaled = function (ratio) {
  return this.getTop(this.height * ratio);
}

Rect.prototype.getBottomScaled = function (ratio) {
  return this.getBottom(this.height * ratio);
}

Rect.prototype.getLeftScaled = function (ratio) {
  return this.getLeft(this.width * ratio);
}

Rect.prototype.getRightScaled = function (ratio) {
  return this.getRight(this.width * ratio);
}

// --

Rect.prototype.getTopPart = function (count, index) {
  return new Rect(
    this.left, this.top + (this.height / count) * index,
    this.width, this.height / count);
}

Rect.prototype.getLeftPart = function (count, index) {
  return new Rect(
    this.left + (this.width / count) * index, this.top,
    this.width / count, this.height);
}

Rect.prototype.getGridPart = function (xCount, yCount, x, y) {
  return new Rect(
    this.left + (this.width / xCount) * x, this.top + (this.height / yCount) * y,
    this.width / xCount, this.height / y);
}

Rect.prototype.inflate = function (dx, dy) {
  return new Rect(
    this.left - dx, this.top - dy,
    this.width + 2 * dx, this.height + 2 * dy);
}

Rect.prototype.inflateSize = function (dw, dh) {
  return new Rect(this.left, this.top, this.width + dw, this.height + dh);
}

Rect.prototype.inflateRel = function (dxRatio, dyRatio) {
  return this.inflate(this.width * (dxRatio - 1), this.height * (dyRatio - 1));
}

// --

Rect.prototype.fitSquare = function () {
  const size = Math.min(this.width, this.height);
  return new Rect(
    this.left + (this.width - size) / 2, this.top + (this.height - size) / 2, size, size);
}

/* From aparamonov
Rect<T> getBottom(T height) => new Rect(left, bottom - height, width, height);
Rect<T> getTop(T height) => new Rect(left, top, width, height);
Rect<T> getLeft(T width) => new Rect(left, top, width, height);
Rect<T> getRight(T width) => new Rect(right - width, top, width, height);
Rect<T> getTopLeft(T width, T height) => new Rect(left, top, width, height);
Rect<T> getTopRight(T width, T height) => new Rect(right - width, top, width, height);
Rect<T> getBottomRight(T width, T height) => new Rect(right - width, bottom - height, width, height);

Rect<T> cutLeft(T dw) => new Rect(left + dw, top, width - dw, height);
Rect<T> cutTop(T dh) => new Rect(left, top + dh, width, height - dh);
Rect<T> cutBottom(T dh) => new Rect(left, top, width, height - dh);
Rect<T> cutRight(T dw) => new Rect(left, top, width - dw, height);

Rect<T> below(T height) => new Rect(left, bottom, width, height);
Rect<T> above(T height) => new Rect(left, top - height, width, height);

Rect<T> toTheLeft(T width) => new Rect(left - width, top, width, height);
Rect<T> toTheRight(T width) => new Rect(right, top, width, height);
Rect<T> toTheTop(T height) => above(height);
Rect<T> toTheBottom(T height) => below(height);

Rect<T> getTopScaled(num ratio) => getTop(height * ratio);
Rect<T> getBottomScaled(num ratio) => getBottom(height * ratio);
Rect<T> getLeftScaled(num ratio) => getLeft(width * ratio);
Rect<T> getRightScaled(num ratio) => getRight(width * ratio);

Rect<T> getTopPart(int count, int index) =>
    new Rect(left, top + (height / count) * index, width, height / count);
Rect<T> getLeftPart(int count, int index) =>
    new Rect(left + (width / count) * index, top, width / count, height);
Rect<T> getGridPart(int xCount, int yCount, int x, int y) =>
    new Rect(left + (width / xCount) * x, top + (height / yCount) * y,
             width / xCount, height / yCount);

Rect<T> inflate(num dx, num dy) => new Rect(left - dx, top - dy, width + 2 * dx, height + 2 * dy);
Rect<T> inflateSize(num dw, num dh) => new Rect(left, top, width + dw, height + dh);
Rect<T> inflateRel(num dxRatio, num dyRatio) => inflate(width * (dxRatio - 1), height * (dyRatio - 1));
/**/

