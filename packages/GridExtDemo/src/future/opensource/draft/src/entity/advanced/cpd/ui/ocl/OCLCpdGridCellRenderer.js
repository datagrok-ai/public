import {ContainerGridCellRenderer} from "../../../../ui/ContainerGridCellRenderer";
import {OCLCpdEntityRenderer} from "./OCLCpdEntityRenderer";

export class OCLCpdGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new OCLCpdEntityRenderer());
    }
}