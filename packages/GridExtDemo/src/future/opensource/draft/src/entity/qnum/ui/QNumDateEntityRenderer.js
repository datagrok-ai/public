import {AbstractVertLayoutTextRenderer} from "../../ui/AbstractVertLayoutTextRenderer";
import {FZModifier} from "../../../lang/FZModifier";
import {TextUtils} from "../../../utils/TextUtils";
import {DateUtils} from "../../../utils/DateUtils";

export class QNumDateEntityRenderer extends AbstractVertLayoutTextRenderer
{
    fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors) {

        if(entity === null)
            return;

        const nTime = entity.getTime();

        let strDate = "";
        if(!isNaN(nTime))
          strDate = TextUtils.formatDate(nTime);

        const modifierFuzzy = entity.getFuzzyModifier();
        const strModifier = modifierFuzzy === null ? "" : (modifierFuzzy === FZModifier.EXACT_EQUAL_TO ? "" : modifierFuzzy.toString());
        const strMainValue =  entity.getValue() === null || entity.getValue() === DG.FLOAT_NULL ? "" : strModifier + TextUtils.formatNumber(entity.getValue());

        arTextLabels.push(strMainValue);
        arTextLabels.push(strDate);

        const font = this.getFont();
        const nHFont = TextUtils.getFontSize(font);
        const nHFontSub = nHFont < 0 ? nHFont : nHFont -2;
        const fontSub = TextUtils.setFontSize(font, nHFontSub);

        arTextFonts.push(font);
        arTextFonts.push(fontSub);

        arTextColors.push("black");
        arTextColors.push("black");

        arBackColors.push(null);
        arBackColors.push(null);
    }

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        const nTime = entity.getTime();
        if(DateUtils.isRecent(nTime, DateUtils.LATER_DAY_COUNT))
        {
            const nSide = Math.min(7, Math.min(nW,nH));
            g.fillStyle = 'green';
            g.fillRect(nX + nW - nSide, nY, nSide, nSide);
        }
    }
}
