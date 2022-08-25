import {AbstractCpdEntityRenderer} from "../AbstractCpdEntityRenderer";
import {createRDKIt} from "../../../../../chem/rdkit/RDKit_minimal_2021.03_02";

export class RDKitCpdEntityRenderer extends AbstractCpdEntityRenderer
{
     constructor(rdKitModule)
     {
        super();
        this.m_rdKitModule = rdKitModule;
     }

    paintStructure(g, entity, nX, nY, nW, nH, crBack) {

        const strSmiles = entity.getSmiles();
        const mol = this.m_rdKitModule.get_mol(strSmiles);

        g.save();

        const opts = {
            "clearBackground": false,
            "offsetx": nX, "offsety": -nY,
            "width": nW,
            "height": nH,
            "bondLineWidth": 1,
            "minFontSize": 9,
            "highlightBondWidthMultiplier": 12,
            "dummyIsotopeLabels": false,
            "atomColourPalette": {
                16: [0.498, 0.247, 0.0],
                9: [0.0, 0.498, 0.0],
                17: [0.0, 0.498, 0.0],
            }
        };
          /*
        let opts = {
            "clearBackground": false,
            "offsetx": nX, "offsety": nY,
            "width": Math.floor(nW),
            "height": Math.floor(nH),
            "bondLineWidth": 1,
            //"fixedScale": 0.07,
            "minFontSize": 9,
            "highlightBondWidthMultiplier": 12,
            "dummyIsotopeLabels": false,
            "atomColourPalette": {
                16: [0.498, 0.247, 0.0],
                9: [0.0, 0.498, 0.0],
                17: [0.0, 0.498, 0.0],
            },
        };  */


        mol.draw_to_canvas_with_highlights(g.canvas, JSON.stringify(opts));
        //mol.draw_to_canvas_with_offset(g.canvas, nX, -nY, nW, nH);

        mol.delete();
        g.restore();
    }
}


RDKitCpdEntityRenderer.create = async function()
{
    const rdKitModule = await createRDKIt("6976ad1747eb075e7a1d.module.wasm");

    const renderer = new RDKitCpdEntityRenderer(rdKitModule);
    return renderer;
}
