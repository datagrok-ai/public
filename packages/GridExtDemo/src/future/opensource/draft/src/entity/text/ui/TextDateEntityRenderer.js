import {AbstractVertLayoutTextRenderer} from "../../ui/AbstractVertLayoutTextRenderer";
import {DateUtils} from "../../../utils/DateUtils";
import {TextUtils} from "../../../utils/TextUtils";

export class TextDateEntityRenderer extends AbstractVertLayoutTextRenderer
{
    fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors) {

        if(entity === null)
            return;

        const nTime = entity.getTime();

        let strDate = "";
        if(!isNaN(nTime))
            strDate = TextUtils.formatDate(nTime);

        const strText = entity.getText();

        arTextLabels.push(strText);
        arTextLabels.push(strDate);

        const font = this.getFont();
        const nHFont = TextUtils.getFontSize(font);
        const nHFontSub = nHFont < 0 ? nHFont : nHFont -2;
        const fontSub = TextUtils.setFontSize(font, nHFontSub);

        arTextFonts.push(font);
        arTextFonts.push(fontSub);

        arTextColors.push("black");
        arTextColors.push("black");

        arBackColors.push("white");
        arBackColors.push("white");
    }

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        const nTime = entity.getTime();
        if(DateUtils.isRecent(entity, DateUtils.LATER_DAY_COUNT))
        {
            const nSide = Math.min(7, Math.min(nW,nH));
            g.fillStyle = 'green';
            g.fillRect(nX + nW - nSide, nY, nSide, nSide);
        }
    }
}
