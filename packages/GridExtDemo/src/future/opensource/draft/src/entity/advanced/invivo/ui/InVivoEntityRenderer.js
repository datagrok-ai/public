import {SemEntityRenderer} from "../../../ui/SemEntityRenderer";
import {TextUtils} from "../../../../utils/TextUtils";

export class InVivoEntityRenderer extends SemEntityRenderer {

    getPreferredCellWidth() {return 50;}

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        let mapInvivo = entity;
        if(!mapInvivo.hasValue())
            return;

        const font = this.getFont();
        g.font = font;//"12px Arial";
        g.fillStyle = "black";
        g.strokeStyle = "black";
        g.textAlign = "center";
        g.textBaseline = "top";

        let str = "InVivo";
        str = TextUtils.trimText(str, g, nW);
        let tm = g.measureText(str);
        let nHText = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent;
        let nYShift = (nH - nHText)/2;

        g.fillText(str, nX + nW / 2, nY + nYShift);
    }
}
