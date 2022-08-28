export class MathUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}


MathUtils.isNullValue = function (obValue)
{
    if(typeof obValue === "number" && isNaN(obValue))
        return true;

    return obValue === undefined || obValue === null || obValue === "" || obValue === DG.FLOAT_NULL || obValue === DG.INT_NULL;
}


MathUtils.generateSubsequentIntegerArray = function(nValFr, nLength)
{
 const ar = new Array(nLength);
 for(var n=0; n<nLength; ++n)
 {
     ar[n] = nValFr + n;
 }

 return ar;
}