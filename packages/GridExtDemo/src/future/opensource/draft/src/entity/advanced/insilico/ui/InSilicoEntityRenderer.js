import {KeyValueRenderer} from "../../../ui/KeyValueRenderer";

export class InSilicoEntityRenderer extends KeyValueRenderer {

    getPreferredCellWidth()
    {
        const i = this.getInsets();
        return 90 + i.getL() + i.getR();
    }

    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 90 + i.getT() + i.getB();
    }


    preparePairsAndUI(entity, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColors,
                      arValsFonts, arValsForeColors, arValsBackColors)
    {
        let map = entity;
        if (map === null || map === undefined)
            return;


        const font = this.getFont();//"13px Dialog";
        const crBad = "pink";
        const crGood = "lightgreen";
        const crOther = "yellow";

        arKeys.push("Lipinski Violations");
        arKeys.push("MW");
        arKeys.push("HBD");
        arKeys.push("HBA");
        arKeys.push("CLogP");
        arKeys.push("PSA");
        arKeys.push("NRotB");


        let str = map.getLipinskiViolations().toString();
        arVals.push(str);

        str = map.getMW().toString();
        arVals.push(str);

        str = map.getHBD().toString();
        arVals.push(str);

        str = map.getHBA().toString();
        arVals.push(str);

        str = map.getCLOGP().toString();
        arVals.push(str);

        str = map.getTPSA().toString();
        arVals.push(str);

        str = map.getNROTB().toString();
        arVals.push(str);

        for(var nKey=0; nKey<arKeys.length; ++nKey)
        {
            arKeysFonts.push(font);
            arValsFonts.push(font);
            arKeysForeColors.push("black");
            arValsForeColors.push("black");
            arKeysBackColors.push("lightGray");
        }


        let cr = map.getLipinskiViolations() === -1 ? "white" : map.getLipinskiViolations() > 2 ? crBad : map.getLipinskiViolations() > 0 ? crOther : crGood;
        arValsBackColors.push(cr);

        cr = isNaN(map.getMW()) ? "white" : (map.getMW() > 500 || map.getMW() < 150) ? crBad : crGood;
        arValsBackColors.push(cr);

        cr = map.getHBD() === -1 ? "white" : map.getHBD() > 5 ? crBad : crGood;
        arValsBackColors.push(cr);

        cr = map.getHBA() === -1 ? "white" : map.getHBA() > 10 ? crBad : crGood;
        arValsBackColors.push(cr);

        cr = isNaN(map.getCLOGP()) ? "white" : (map.getCLOGP() > 5 || map.getCLOGP() < -0.7) ? crBad : crGood;
        arValsBackColors.push(cr);

        cr = map.getTPSA() === -1 ? "white" : map.getTPSA() > 150 ? crBad : crGood;
        arValsBackColors.push(cr);

        cr = isNaN(map.getNROTB()) ? Color.white : map.getNROTB() > 10 ? crBad : crGood;
        arValsBackColors.push(cr);
    }
}
