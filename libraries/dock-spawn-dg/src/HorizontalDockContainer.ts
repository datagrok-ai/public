import { SplitterDockContainer } from "./SplitterDockContainer";
import { Utils } from "./Utils";
import { DockManager } from "./DockManager";
import { ContainerType } from "./ContainerType";
import { IDockContainer } from "./interfaces/IDockContainer";

export class HorizontalDockContainer extends SplitterDockContainer {

    constructor(dockManager: DockManager, childContainers: IDockContainer[]) {
        super(Utils.getNextId('horizontal_splitter_'), dockManager, childContainers, false)
        this.containerType = ContainerType.horizontal;
    }
}
