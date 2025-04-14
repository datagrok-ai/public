import { ISize } from "./ISize";
import { IDockContainer } from "./IDockContainer";

export interface IDockContainerWithSize extends IDockContainer {
    state: ISize;
}
