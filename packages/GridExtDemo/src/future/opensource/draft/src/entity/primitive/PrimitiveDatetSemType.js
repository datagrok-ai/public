import {PrimitiveSemType} from "./PrimitiveSemType";

export class PrimitiveDatetSemType extends PrimitiveSemType
{
    constructor()
    {
        super("images/date.png");
    }
}

PrimitiveDatetSemType.Instance = new PrimitiveDatetSemType();
