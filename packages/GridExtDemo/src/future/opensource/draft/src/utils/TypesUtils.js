export class TypesUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}

TypesUtils.isString = function(ob)
{
    const b = TypesUtils.isTypeOf(ob, String);
    return b;
}

TypesUtils.isNumber = function(ob)
{
    const b = TypesUtils.isTypeOf(ob, Number);
    return b;
}


TypesUtils.isTypeOf = function(ob, obType)
{
    const b = typeof ob === obType.name.toLowerCase();
    return b;
}
