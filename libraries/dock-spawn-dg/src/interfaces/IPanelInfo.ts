import { ContainerType } from "../ContainerType";
import { Point } from "../Point";
import { IState } from "./IState";

export interface IPanelInfo {
    containerType: ContainerType;
    state: IState;
    position: Point;
    isHidden: boolean;
}
