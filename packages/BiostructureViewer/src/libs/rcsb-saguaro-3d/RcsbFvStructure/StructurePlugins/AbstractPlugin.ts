import {RcsbFvSelectorManager} from "../../RcsbFvSelection/RcsbFvSelectorManager";

export class AbstractPlugin {
    protected readonly selection: RcsbFvSelectorManager;

    constructor(selection: RcsbFvSelectorManager) {
        this.selection = selection;
    }
}
