/// this file was generated automatically from d4 classes declarations
import { toDart } from "../wrappers";
let api = <any>window;

export function renderMultipleHistograms(g: CanvasRenderingContext2D, bounds: any, histograms: Array<Int32List>, options?: {categoryColumn?: any, colors?: Array<number>, tension?: number, normalize?: boolean, markerSize?: number, fill?: boolean, minBin?: number, maxBin?: number, localMaximum?: boolean}): any
  { return api.grok_renderMultipleHistograms(toDart(g), toDart(bounds), toDart(histograms), toDart(options?.categoryColumn), toDart(options?.colors), toDart(options?.tension), toDart(options?.normalize), toDart(options?.markerSize), toDart(options?.fill), toDart(options?.minBin), toDart(options?.maxBin), toDart(options?.localMaximum)); }

