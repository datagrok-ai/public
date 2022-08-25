import {PrimitiveSemType} from "./PrimitiveSemType";

export class PrimitiveStringSemType extends PrimitiveSemType
{
    constructor()
    {
        super("images/string.png");
    }
}

PrimitiveStringSemType.Instance = new PrimitiveStringSemType();
