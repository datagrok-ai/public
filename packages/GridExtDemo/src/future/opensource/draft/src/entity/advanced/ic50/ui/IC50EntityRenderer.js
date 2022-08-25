import {AbstractVertLayoutTextRenderer} from "../../../ui/AbstractVertLayoutTextRenderer";
import {FZModifier} from "../../../../lang/FZModifier";
import {TextUtils} from "../../../../utils/TextUtils";
import {DateUtils} from "../../../../utils/DateUtils";

export class IC50EntityRenderer extends AbstractVertLayoutTextRenderer
{
    fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors)
    {
        if(entity === null)
        {
         return;
        }

        const nTime = entity.getLastTested();

        let strDate = null;
        let crLastTested = null;
        if (nTime !== Number.NaN){

            strDate = TextUtils.formatDate(nTime);
            crLastTested = DateUtils.getRecentDateColor(nTime);//DateUtils.isRecent(nTime, DateUtils.LATEST_DAY_COUNT) ? DateUtils.LATEST_DAY_COLOR : DateUtils.isRecent(nTime, DateUtils.LATER_DAY_COUNT) ? DateUtils.LATER_DAY_COLOR : null;
        }

        const strFuzzy = entity.getFuzzyModifier() === FZModifier.EXACT_EQUAL_TO ? "" : entity.getFuzzyModifier().toString();

        arTextLabels.push(strFuzzy + entity.getMedianIC50().toString() + "±" + entity.getStdv().toString());
        arTextLabels.push(entity.getPctEff());
        arTextLabels.push(strDate);

        const font = this.getFont();
        const nHFont = TextUtils.getFontSize(font);
        const nHFontSub = nHFont < 0 ? nHFont : nHFont -1;
        const fontSub = TextUtils.setFontSize(font, nHFontSub);

        arTextFonts.push(font);
        arTextFonts.push(fontSub);
        arTextFonts.push(fontSub);

        arTextColors.push("black");
        arTextColors.push("black");
        arTextColors.push("black");

        arBackColors.push("white");
        arBackColors.push(IC50EntityRenderer.m_crPctEff);
        arBackColors.push(crLastTested);
    }

    paint(g, entity, nX, nY, nW, nH, crBack) {

        const fVal = entity.getMedianIC50();
        if(fVal !== undefined) {
            if (fVal > 10.0)
                crBack = "LightPink";
            else if (fVal < 2.0)
                crBack = "LightGreen";
            else
                crBack = "Khaki";
        }
        super.paint(g, entity, nX, nY, nW, nH, crBack);

        if(entity === null)
        {
            return;
        }


        const nTime = entity.getLastTested();
        if(DateUtils.isRecent(nTime, DateUtils.LATER_DAY_COUNT))
        {
            const nSide = Math.min(7, Math.min(nW,nH));
            g.fillStyle = 'green';
            g.fillRect(nX + nW - nSide, nY, nSide, nSide);
        }
    }
}


IC50EntityRenderer.m_crPctEff = "PaleTurquoise";
IC50EntityRenderer.m_crLatest = "LightPink";
IC50EntityRenderer.m_crLater = "cyan";
