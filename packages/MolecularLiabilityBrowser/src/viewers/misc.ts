import * as DG from 'datagrok-api/dg';

export class MiscMethods {
  static NoSchemeItem: string = '<not selected>';

  /** processes JSON to derive scheme names */
  static extractSchemes(json: any) {
    const rawSchemeNames = Object.keys(json.cdr_ranges);
    const schemesLst = [MiscMethods.NoSchemeItem];
    rawSchemeNames.forEach((str) => {
      const strArr = str.split('_');
      if (schemesLst.includes(strArr[0]) === false)
        schemesLst.push(strArr[0]);
    });
    return schemesLst;
  }

  // color interpolation
  static interpolateColors(color1: string, color2: string, steps: number) {
    function interpolateColor(color1: number[], color2: number[], factor: number) {
      if (arguments.length < 3)
        factor = 0.5;


      const result = color1.slice();
      for (let i = 0; i < 3; i++)
        result[i] = Math.round(result[i] + factor * (color2[i] - color1[i]));


      const hexCol = '#' + ((1 << 24) + (result[0] << 16) + (result[1] << 8) + result[2]).toString(16).slice(1);
      return hexCol;
    };

    const stepFactor = 1 / (steps - 1);
    const interpolatedColorArray = [];
    const colors1 = color1.match(/\d+/g)!.map(Number);
    const colors2 = color2.match(/\d+/g)!.map(Number);

    for (let i = 0; i < steps; i++)
      interpolatedColorArray.push(interpolateColor(colors1, colors2, stepFactor * i));

    return interpolatedColorArray;
  }

  // ---- Resizing ----
  static setDockSize(view: DG.TableView, node: DG.DockNode, nodeContent: HTMLElement | DG.Viewer): DG.DockNode {
    const rootNodeHeight: number = view.dockManager.rootNode.container.containerElement.clientHeight;

    //@ts-ignore
    const newHeight = $('#feature-viewer').outerHeight(true) + 57;
    const ratioHeight: number = rootNodeHeight != 0 ? 1 / (rootNodeHeight / newHeight) : 0.5;

    //@ts-ignore
    return view.dockManager.dock(nodeContent, 'down', node, 'Sequence', ratioHeight);
  }
}
