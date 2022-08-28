export class RandUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}
   /*
RandUtils.generateSubsequentIntegerArray = function(nValFr, nLength)
{
    const ar = new Array(nLength);
    for(var n=0; n<nLength; ++n)
    {
        ar[n] = nValFr + n;
    }

    return ar;
}    */


RandUtils.generateRandomInteger = function(nIntFr, nIntTo) {

    let nInt = Math.floor(Math.random() * (nIntTo - nIntFr + 1)) + nIntFr;
    return nInt;
}

RandUtils.generateRandomDate = function(nYearFr, nYearTo) {

    if(nYearFr > nYearTo)
        throw new Error("Starting year (" + nYearFr + ") cannot exceed ending year (" + nYearTo + ")");

    const nD =  Math.floor(Math.random() * 28);
    const nM = Math.floor(Math.random() * 12);
    const nY = Math.floor(Math.random() * (nYearTo - nYearFr + 1)) + nYearFr;

    const nTimeNow = new Date().getTime();
    let nTime = new Date(nY, nM, nD).getTime();
    if(nTime > nTimeNow)
        nTime = nTimeNow;

    return nTime;
}

RandUtils.generateRandomRecentDate = function() {
    const nMonth = 10;//Math.floor(Math.random() * 7);
    const nDay =  1 + Math.floor(Math.random() * 8);
    const nTimeLastTested = new Date(2021, nMonth, nDay).getTime();
    return nTimeLastTested;
}

class SomeObject {
    constructor() {
        this.m_value = Math.random();
    }
}


RandUtils.generateRandomObjectDataFrame = function(nRecordCount) {
    let nColCount = 50;
    let arCols = new Array(nColCount);

    var nCol=0;
    for(; nCol<nColCount; ++nCol)
    {
        arCols[nCol] = DG.Column.fromType(DG.TYPE.OBJECT, "Column " + nCol.toString(), nRecordCount);
        // arCols[nCol].semType = IC50SemType.INSTANCE;//IC50Map.semType;
    }


    for(var n=0; n<nRecordCount; ++n)
    {

        for(nCol=0; nCol<nColCount; ++nCol)
        {
            arCols[nCol].set(n, new SomeObject());
        }
    }

    let dframe =  null;
    try{dframe = DG.DataFrame.fromColumns(arCols);}
    catch(e)
    {
        let a = 0;
    }

    return dframe;
}


RandUtils.generateRandomCpdDataFrame = function(nRecordCount, obPackage) {
    if(nRecordCount === undefined)
        throw new Error("The number of records cannot be undefined.");


    let nColCount = 2;
    let arCols = new Array(nColCount);
    arCols[0] = DG.Column.fromType(DG.TYPE.OBJECT, "Cpd", nRecordCount);
    arCols[0].semType = "Cpd";

    arCols[1] = DG.Column.fromType(DG.TYPE.STRING, "Mol", nRecordCount);
    arCols[1].semType = "Molecule";

    let cpd = null;
    for(var n=0; n<nRecordCount; ++n) {
        try {
            cpd = new CpdEntity();
            arCols[0].set(n, cpd);
            arCols[1].set(n, cpd.getSmiles());
        } catch (e) {
            let rdfg = 0;
        }
    }

    let dframe =  null;
    try{dframe = DG.DataFrame.fromColumns(arCols);}
    catch(e)
    {
        let a = 0;
    }

    return dframe;
}


