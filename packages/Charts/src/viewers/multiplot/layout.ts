export function MPlotLayout2(visibleIndex: number, viewerThis: any) {
  return 0;
}

export class MPLayout {
  constructor() {}

  // Assign size and position of viewer's parts (charts, controls)
  // input: height of viewer's container
  // also takes data from this.plots:
  //    [{height: '10px', title: 'title0'},{height: '20%'}, {height: 'free'}]
  // output: [{top: '20px', height: '10px', topTitle: '0px', titleHeight}, ...]
  parseHeight(heightAll: number, plots : any[], layoutSettings : any): any[] {
    const r: any[] = [];
    const height = heightAll - layoutSettings.globalMarginBottom -
        layoutSettings.globalMarginTop;
    let floatHeight = height;
    let totalFlexPoints = 0;

    for (let i = 0; i < plots.length; i++) { // first cycle
      const titleHeight = plots[i].title ?
        parseFloat(layoutSettings.defaultTitleHeight) : 0;
      let plotHeight = 0;
      let flexPoints = 0;
      if (plots[i].height.includes('px'))
        plotHeight = parseFloat(plots[i].height);

      if (plots[i].height.includes('%'))
        plotHeight = parseFloat(plots[i].height) * height / 100;

      if (plots[i].height.includes('flex') && plots[i].show) {
        flexPoints = parseFloat(plots[i].height);
        totalFlexPoints += flexPoints;
      };
      if (plots[i].show == 0) plotHeight = 0;

      r.push({
        flexPoints: flexPoints,
        height: plotHeight,
        titleHeight: titleHeight,
        titleText: plots[i].title,
      });
      floatHeight -= plotHeight + titleHeight +
        1 * (layoutSettings.plotTitleHighMargin + layoutSettings.plotTitleLowMargin);
    } // first cycle

    // weight (in pixels) of one flex point
    const flexPointHeight = floatHeight / totalFlexPoints;

    let currentTop = layoutSettings.globalMarginTop;
    let echartIndex = 0; // index in echarts library
    // second cycle (set positions and height)
    for (let j = 0; j < r.length; j++) {
      currentTop += layoutSettings.plotTitleHighMargin;
      r[j].titleTop = currentTop;
      currentTop += r[j].titleHeight;
      if (plots[j].show) {
        currentTop += layoutSettings.plotTitleLowMargin;
        r[j].top = currentTop;
        if (r[j].flexPoints > 0.00001)
          r[j].height = r[j].flexPoints * flexPointHeight;
        ;
        currentTop += r[j].height;
        plots[j].i = echartIndex;
        plots[j].floatHeight = r[j].height;
        echartIndex++;
      };
    }

    // correct height of the last visible plot to avoid overlap of controls
    let j = plots.length - 1;
    while (j > -1 && !plots[j].show)
      j--;

    if (j < plots.length - 1 && j >= 0)
      r[j].height += layoutSettings.headerHiddenShift;

    return r;
  } // parseHeight
}

