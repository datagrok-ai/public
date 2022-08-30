import {PrimitiveSemType} from "./PrimitiveSemType";

export class PrimitiveQNumSemType extends PrimitiveSemType
{
    constructor()
    {
        super("images/qnum.png");
    }
}

PrimitiveQNumSemType.Instance = new PrimitiveQNumSemType();
