export type ChemspaceResult = {
    csId: string;
    link: string;
    smiles?: string;
    molFormula?: string;
    cas?: string;
    mfcd?: string;
    properties: ChemspaceProperties;
    offerCount: number;
    offers: ChemspaceOffer[];
}

export type ChemspaceProperties = {
    mw: number;
    hac: number;
    logp: number;
    rotb: number;
    hba: number;
    hbd: number;
    ringCount: number;
    fsp3: number;
    tpsa: number;
}

export type ChemspaceOffer = {
    vendorName: string;
    vendorCode: string;
    leadTimeDays?: number;
    purity: number;
    prices: ChemspacePrice[];
}

export type ChemspacePrice = {
    packMg: number;
    priceUsd?: number;
    priceEur?: number;
}

export type ChemspacePricesTableItem = {
    packMg: number;
    priceUsd?: number;
    priceEur?: number;
    vendorName: string;
    vendorCode: string;
    leadTimeDays?: number;
    purity: number;
}