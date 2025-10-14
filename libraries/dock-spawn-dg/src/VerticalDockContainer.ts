import { SplitterDockContainer } from "./SplitterDockContainer";
import { Utils } from "./Utils";
import { ContainerType } from "./ContainerType";
import { DockManager } from "./DockManager";
import { IDockContainer } from "./interfaces/IDockContainer";

export class VerticalDockContainer extends SplitterDockContainer {

    constructor(dockManager: DockManager, childContainers: IDockContainer[]) {
        super(Utils.getNextId('vertical_splitter_'), dockManager, childContainers, true)
        this.containerType = ContainerType.vertical;
    }
}
