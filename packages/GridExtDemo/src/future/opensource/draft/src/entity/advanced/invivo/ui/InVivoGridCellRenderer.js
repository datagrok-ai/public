import {ContainerGridCellRenderer} from "../../../ui/ContainerGridCellRenderer";
import {InVivoEntityRenderer} from "./InVivoEntityRenderer";

export class InVivoGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new InVivoEntityRenderer());
    }
}