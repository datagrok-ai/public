export class TextDimensionsCalculator {
  private static instance: TextDimensionsCalculator;
  private canvas: HTMLCanvasElement;
  private constructor() { }

  private static getInstance(): TextDimensionsCalculator {
    if (!TextDimensionsCalculator.instance) {
      TextDimensionsCalculator.instance = new TextDimensionsCalculator();
      TextDimensionsCalculator.instance.canvas = document.createElement('canvas');
    }
    return TextDimensionsCalculator.instance;
  }

  static getTextDimensions(text: string, fontSize: number): {width: number, height: number} {
    const canvas = TextDimensionsCalculator.getInstance().canvas;
    const context = canvas.getContext('2d');
    if (!context)
      throw new Error('Canvas 2D context is not available');

    context.font = `${fontSize}px Arial`;
    const scaleFactor = 1.1;

    return {
      width: context.measureText(text).width * scaleFactor,
      height: fontSize * scaleFactor
    };
  }
}

