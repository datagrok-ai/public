import {SemType} from "./SemType";

export class NullSemType extends SemType
{
    constructor() {
        super(undefined, []);
    }
}