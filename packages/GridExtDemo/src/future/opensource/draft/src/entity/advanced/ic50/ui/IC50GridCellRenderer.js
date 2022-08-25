import {IC50EntityRenderer} from "./IC50EntityRenderer";
import {ContainerGridCellRenderer} from "../../../ui/ContainerGridCellRenderer";

export class IC50GridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new IC50EntityRenderer());
    }
}