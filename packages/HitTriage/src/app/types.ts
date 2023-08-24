import * as DG from 'datagrok-api/dg';

export type IDescriptorTree = {
    [key: string]: {
      descriptors: Array<Descriptor>,
    } & Descriptor;
}

export type Descriptor = {
    name: string,
    description: string,
};

export type IFunctionArgs<T = any> = {
    [key: string]: T,
}

export type HitTriageTemplateFunction = {
    package: string,
    name: string,
    args: IFunctionArgs,
}

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: IFunctionArgs
    }
}

export type HitTriageTemplate = {
    name: string,
    key: string,
    compute: HitTriageTemplateCompute,
    submit?: HitTriageTemplateSubmit,
    campaignFields: HitTriageCampaignField[],
    dataSourceType: IngestType,
}

export const CampaignFieldTypes = {
  String: DG.TYPE.STRING,
  Number: DG.TYPE.FLOAT,
  Boolean: DG.TYPE.BOOL,
  Date: DG.TYPE.DATE_TIME,
} as const;

export type HitTriageCampaignFieldType = keyof typeof CampaignFieldTypes;

export type HitTriageCampaignField = {name: string, type: HitTriageCampaignFieldType, required?: boolean};

export type IngestType = 'File' | 'Query';

export type HitTriageTemplateIngest = {
    type: IngestType,
    query: string,
    molColName: string,
};

export type HitTriageTemplateCompute = {
    descriptors: {
        enabled: boolean,
        args: string[],
    }
    functions: HitTriageTemplateFunction[],
};

export type HitTriageTemplateSubmit = {
    fName: string,
    package: string
};

export type HitTriageCampaignStatus = 'In Progress' | 'Submitted';

export type HitTriageCampaign = {
    name: string,
    templateName: string,
    status: HitTriageCampaignStatus,
    createDate: string,
    campaignFields: {[key: string]: any},
    filters: {[key: string]: any}[],
    ingest: HitTriageTemplateIngest,
};

export type IChemFunctionsDialogResult = {
    okProxy: () => void,
    root: HTMLElement,
};

export type INewTemplateResult<T> = {
    template: Promise<T>,
    root: HTMLElement,
    cancelPromise: Promise<void>,
}

// ##################### HIT DESIGN TYPES #####################

export type HitDesignTemplate = Omit<HitTriageTemplate, 'dataSourceType'> & {stages: string[]};

export type HitDesignCampaign = Omit<HitTriageCampaign, 'filters' | 'ingest'>;
