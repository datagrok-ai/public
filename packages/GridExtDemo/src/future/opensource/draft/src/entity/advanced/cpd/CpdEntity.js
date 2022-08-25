import {SemEntity} from "../../SemEntity";
import {RandUtils} from "../../../utils/RandUtils";

export class CpdEntity  extends SemEntity
{
    constructor(strSmiles, bMolFile, strSampleId, strConceptId, nTimeLastUpdated, strChemist)
    {
        super();

        if(strSmiles === undefined)
            strSmiles = CpdEntity.SAMPLE_SMILES[CpdEntity.INSTAMCE%CpdEntity.SAMPLE_SMILES.length];

        this.m_bMolFile = false;
        if(bMolFile!== undefined)
            this.m_bMolFile = bMolFile;

        let strGnfId = "GNF-XXX-XXXX-" + RandUtils.generateRandomInteger(1, 100);
        let nCpdSidId = RandUtils.generateRandomInteger(1000008, 9653348);

        if(strSampleId === undefined)
            strSampleId = "NVP-ASF-" + RandUtils.generateRandomInteger(10000, 90000);

        this.m_strSmiles = strSmiles;
        this.m_strGnfRegId  = strGnfId;

        this.m_strSampleId = strSampleId;
        this.m_strConceptId = strConceptId;
        this. m_nCpdSidId = nCpdSidId;

        this.m_strLegacyId = null;
        this.m_strSrcCodeId = null;

        if(nTimeLastUpdated === undefined)
            nTimeLastUpdated = NaN;//new Date().getTime();//RandUtils.generateRandomDate(2007, 2021);

        this.m_nLastTested = nTimeLastUpdated;

        if(strChemist === undefined)
            strChemist = "chemist";

        this.m_strChemist = strChemist;


        ++CpdEntity.INSTAMCE;
    }

    getDefaultValue()
    {
        return this.getConceptId();
    }

    isMolFile() {return this.m_bMolFile;}
    getSmiles() {return this.m_strSmiles};
    getGnfRegId() {return this.m_strGnfRegId;}
    getSampleId() {return this.m_strSampleId;}
    getConceptId() {return this.m_strConceptId;}
    getCpdSidId() {return this.m_nCpdSidId;}
    getLegacyId() {return this.m_strLegacyId;}
    getSrcCodeId() {return this.m_strSrcCodeId;}
    getLastTested() {return this.m_nLastTested;}
    getChemist() {return this.m_strChemist;}

    toString() {return this.m_strConceptId;}
}

CpdEntity.INSTAMCE = 0;
CpdEntity.SAMPLE_SMILES =[
    "CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3",
    "COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC",
    "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3",
    "CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2",
    "FC(F)(F)c1ccc(OC2CCNCC2)cc1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCCc3ccccc3",
    "COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO[N+](=O)[O-]",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO",
    "CN1CCC(CC1)Oc2ccc(cc2)C(F)(F)F",
    "COc1cc(C)ccc1OC(=O)C(C)c2ccc(CC(C)C)cc2",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)OCCCc2cccnc2",
    "COc1ccc(\C=N\NC(=N)N)c(Cl)c1OC",
    "Nc1ncnc2c1c(Br)cn2[C@@H]3OC[C@@H](O)[C@H]3O",
    "CNc1ncnc2c1c(I)cn2[C@@H]3O[C@H](C)[C@@H](O)[C@H]3O",
    "CN1CCC(O)(CC1)c2ccccc2",
    "OC(COc1ccccc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
    "OC(COc1ccc(Cl)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
    "OC(COc1ccc(Br)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
    "COC(=O)c1ccccc1OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
    "CCC1(CC)CC(CCNC(=O)c2ccc(OC)cc2)OC1=O"];

CpdEntity.generateRandomSmiles = function()
{
    const strSmiles = CpdEntity.SAMPLE_SMILES[CpdEntity.INSTAMCE%CpdEntity.SAMPLE_SMILES.length];
    ++CpdEntity.INSTAMCE;
    return strSmiles;
}
