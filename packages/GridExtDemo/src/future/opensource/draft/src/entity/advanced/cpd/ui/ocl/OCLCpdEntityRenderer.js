import {initOCLModule} from "../../../../../chem/ocl/openchemlib-full.pretty.js";
import {AbstractCpdEntityRenderer} from "../AbstractCpdEntityRenderer";

export class OCLCpdEntityRenderer extends AbstractCpdEntityRenderer
{
    paintStructure(g, cpd, nX, nY, nW, nH, crBack) {

        let strSmiles = cpd.getSmiles();
        let mol = cpd.m_bMolFile ? OCL.Molecule.fromMolfile(strSmiles) : OCL.Molecule.fromSmiles(strSmiles);
        g.save();
        OCL.StructureView.drawMolecule(g.canvas, mol, undefined, nX, nY, nW, nH, g);
        g.restore();
    }
}