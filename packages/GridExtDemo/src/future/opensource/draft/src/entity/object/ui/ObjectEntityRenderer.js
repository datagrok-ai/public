import {AbstractVertLayoutTextRenderer} from "../../ui/AbstractVertLayoutTextRenderer";
import {TextUtils} from "../../../utils/TextUtils";

export class ObjectEntityRenderer extends AbstractVertLayoutTextRenderer
{
    fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors) {

        if(entity === null)
            return;

        const font = this.getFont();
        arTextFonts.push(font);//13px Dialog");
        arTextColors.push("black");
        arBackColors.push("white");

        const obValue = entity.getDefaultValue();
        if(obValue === undefined || obValue === null)
            return; "";

        let strText = "";
        const fValue = Number(obValue);
        if(!isNaN(fValue))
        {
            strText = TextUtils.formatNumber(fValue)
        }
        else strText =  obValue.toString();

        arTextLabels.push(strText);
    }
}


