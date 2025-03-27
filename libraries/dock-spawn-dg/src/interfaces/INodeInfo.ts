import { ContainerType } from "../ContainerType";
import { IState } from "./IState";

export interface INodeInfo {
    containerType: ContainerType;
    state: IState;
    children: INodeInfo[];
}
