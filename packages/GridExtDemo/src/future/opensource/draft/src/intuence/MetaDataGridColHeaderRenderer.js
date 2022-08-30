import {ButtonGridColumnHeaderRenderer} from "../entity/ui/ButtonGridColumnHeaderRenderer";
import {KeyValueRenderer} from "../entity/ui/KeyValueRenderer";
import {AnalysisLoader} from "../service/nx/analysis/AnalysisLoader";
import {CDFCompositePrimitiveType} from "../service/nx/cdf/composite/CDFCompositePrimitiveType";
import {GridUtils} from "../utils/GridUtils";

class MetaDataKeyValueRenderer extends KeyValueRenderer {

    preparePairsAndUI(data, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColots,
                      arValsFonts, arValsForeColors, arValsBackColors) {

        arKeys.push("Name:")
        arVals.push(data.name);
        arKeys.push("Parent:")
        arVals.push(data.parent);
        arKeys.push("Type:")
        arVals.push(data.type);
        arKeys.push("Multiplicity:")
        arVals.push(data.multiplicity);

        const font = this.getFont();
        const nValCount = arVals.length;
        for(let n=0; n<nValCount; ++n)
        {
            arKeysFonts.push("bold " + font);
            arKeysForeColors.push("black");
            arKeysBackColots.push("PapayaWhip");

            arValsFonts.push(font);
            arValsForeColors.push("black");
            arValsBackColors.push("PapayaWhip");
        }
    }
}



export class MetaDataGridColHeaderRenderer extends ButtonGridColumnHeaderRenderer
{
    constructor() {
        super();

     this.m_rendererTootip = new MetaDataKeyValueRenderer();
    }


    createTootipContent(cell)
    {
        const nW = 250;
        let nH = 100;
        let font = cell.style.font;
        if(font === null)
            font = "13px Roboto, Roboto Local";

        this.setFont("bold " + font);

        const col = cell.tableColumn;

        const descriptor = col === null ? null : col.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
        if(descriptor !== null && descriptor !== undefined)
        {
            const typeCDF = descriptor.getPrimitiveType();
            if (typeCDF instanceof CDFCompositePrimitiveType) {

                //Column Name
                let eCanvas = new OffscreenCanvas(nW, nH);
                let g = eCanvas.getContext('2d');
                const strColName = cell.gridColumn.name;
                const nHName = 5 + GridUtils.calcTextrHeight(g, strColName, "bold " + font, nW, 2);
                nH = nHName + 15;

                g.font = "bold " + font;
                const tm = g.measureText("W");
                const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + 4;
                nH += (4*nHFont);

                eCanvas = ui.canvas(nW, nH);
                g = eCanvas.getContext('2d');

                const bPaintBorderOld = this.getPaintBorder();
                const bPaintBackOld = this.getPaintBackground();
                const bPaintLabelOld = this.getPaintLabel();
                this.setPaintBorder(false);
                this.setPaintBackground(false);
                this.setFilterEnabled(false)
                this.setPaintLabel(false);
                g.fillStyle = "wheat";
                g.fillRect(0,0, nW, nH);
                this.paint(g, cell.gridColumn, 0,0, nW, nHName);

                this.setPaintBorder(bPaintBorderOld);
                this.setPaintBackground(bPaintBackOld);
                this.setFilterEnabled(true);
                this.setPaintLabel(bPaintLabelOld);

                //CDF
                const data = {
                    name: typeCDF.getName(),
                    parent: typeCDF.getParentType().getName(),
                    type:typeCDF.getType(),
                    multiplicity: typeCDF.getMultiplicity().toString()
                };

                let fontTmp = cell.style.font;
                if(fontTmp === null)
                    fontTmp = font;

                this.m_rendererTootip.setKeys2ValsWidthRatio(0.3);
                this.m_rendererTootip.setFont(font);
                this.m_rendererTootip.paint(g, data, 0,nHName + 15, nW, nH-nHName-15);
                return eCanvas;
            }
        }

        const bPaintBorderOld = this.getPaintBorder();
        const bPaintBackOld = this.getPaintBackground();
        const bPaintLabelOld = this.getPaintLabel();
        this.setPaintBorder(false);
        this.setPaintBackground(false);
        this.setFilterEnabled(false)
        this.setPaintLabel(false);

        let eCanvas = new OffscreenCanvas(nW, nH);
        let g = eCanvas.getContext('2d');
        nH = 5 + GridUtils.calcTextrHeight(g, cell.gridColumn.name, "bold " + font, nW, 2);

        eCanvas = ui.canvas(nW, nH);
        g = eCanvas.getContext('2d');

        g.fillStyle = "wheat";
        g.fillRect(0,0, nW, nH);
        this.paint(g, cell.gridColumn, 0,0, nW, nH);

        this.setPaintBorder(bPaintBorderOld);
        this.setPaintBackground(bPaintBackOld);
        this.setFilterEnabled(true);
        this.setPaintLabel(bPaintLabelOld);

        return eCanvas;
    }
}
