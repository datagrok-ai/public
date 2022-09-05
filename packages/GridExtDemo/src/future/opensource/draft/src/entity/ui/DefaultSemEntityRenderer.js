import {AbstractVertLayoutTextRenderer} from "./AbstractVertLayoutTextRenderer";

export class DefaultSemEntityRenderer extends AbstractVertLayoutTextRenderer
{
    fillLabelsAndUI(ctx, arTextLabels, arTextFonts, arTextColors, arBackColors) {

        let cell = ctx.getCell();
        let entity = cell.value;
        let obValue = entity.getDefaultValue();
        let strValue = obValue === undefined || obValue === null ? "" : obValue.toString();

        arTextLabels.push(strValue);
        arTextFonts.push("13px Dialog");
        arTextColors.push("black");
        arBackColors.push("white");
    }
}