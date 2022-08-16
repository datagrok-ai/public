import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class ODEview extends DG.ViewBase {

  mainPlotDiv = ui.box();

  constructor(name : string) {
    super({});
    this.name = name;

    let tabControl = ui.tabControl(null, false);

    tabControl.addPane('Survival Chart', () => ui.splitV([
      ui.splitH([
        ui.box( ui.panel([
          ui.inputs([
            ui.choiceInput('AA', 0.1, [0.1, 0.2, 0.3])
          ])
        ]), { style: { maxWidth: '300px' }}),
        this.mainPlotDiv
      ])
    ]));

    this.updatePlot();
  }

  private updatePlot() {
    ui.setUpdateIndicator(this.mainPlotDiv, true);
    grok.functions.call(
      "Simpkpd:rxodeCommandReal", {
      "inputSD": 3
    }).then((result) => {
      ui.setUpdateIndicator(this.mainPlotDiv, false);
      //@ts-ignore
      updateDivInnerHTML(this.mainPlotDiv, ui.image(`data:image/png;base64,${result[ 'plot' ]}`));
    });
  }
}