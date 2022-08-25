import {ContainerGridCellRenderer} from "../../../../ui/ContainerGridCellRenderer";
import {RDKitCpdEntityRenderer} from "./RDKitCpdEntityRenderer";

export class RDKitCpdGridCellRenderer extends ContainerGridCellRenderer
{
    constructor(rendererEntity)
    {
        super(rendererEntity);
    }
}

RDKitCpdGridCellRenderer.create = async function()
{
    const rendererEntity = await RDKitCpdEntityRenderer.create();
    const renderer = new ContainerGridCellRenderer(rendererEntity);

    return renderer;
}