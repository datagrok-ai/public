import {AbstractVertLayoutTextRenderer} from "./AbstractVertLayoutTextRenderer";

export class TestVertLayoutCellRenderer extends AbstractVertLayoutTextRenderer {
  fillLabelsAndUI(entity : any, arTextLabels : string[], arTextFonts: string[], arTextColors : string[], arBackColors: string[]) : void {
    arTextLabels.push('test 1');
    arTextLabels.push('test 2');
    arTextLabels.push('test 3');
    arTextFonts.push('Roboto 12px');
    arTextFonts.push('Roboto 12px');
    arTextFonts.push('Roboto 12px');
    arTextColors.push('black');
    arTextColors.push('black');
    arTextColors.push('black');
    arBackColors.push('white');
    arBackColors.push('white');
    arBackColors.push('yellow');
  }
}
