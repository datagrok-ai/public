import {SemEntity} from "../../SemEntity";

export class InSilicoEntity extends SemEntity {

    constructor() {

        super();

        this.m_nLipinskiViolations = Math.floor(Math.random()*5);
        this.m_fMW = Math.floor(Math.random() * 700); //570;
        this.m_nHBD = Math.floor(Math.random() * 5);//1;
        this.m_nHBA = Math.floor(Math.random() * 10);//8;
        this.m_fCLOGP = Math.floor(Math.random() * 10);//5.35;
        this.m_fTPSA = Math.floor(Math.random() * 100);//77;
        this.m_nNROTB = Math.floor(Math.random() * 7);//3.0;


    }

    getDefaultValue()
    {
        return this.getLipinskiViolations();
    }


    getLipinskiViolations()
    {
        return this.m_nLipinskiViolations;
    }
    getMW() {return this.m_fMW;}
    getHBA() {return this.m_nHBA;}
    getHBD() {return this.m_nHBD;}
    getCLOGP() {return this.m_fCLOGP;}
    getTPSA() {return this.m_fTPSA;}
    getNROTB() {return this.m_nNROTB;}
}
