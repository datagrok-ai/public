import {ContainerGridCellRenderer} from "../../../ui/ContainerGridCellRenderer";
import {InSilicoEntityRenderer} from "./InSilicoEntityRenderer";

export class InSilicoGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new InSilicoEntityRenderer());
    }
}