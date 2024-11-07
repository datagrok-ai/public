import * as DG from 'datagrok-api/dg';

export type AppName = 'Hit Triage' | 'Hit Design' | 'PeptiHit';

export type CampaignsType = {
    'Hit Triage': HitTriageCampaign,
    'Hit Design': HitDesignCampaign,
    'PeptiHit': HitDesignCampaign,
}

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

export type HitTriageTemplateScript = {
    name: string,
    args: IFunctionArgs,
    id: string,
}

export type HitTriageComputeQuery = HitTriageTemplateScript & {
    inputName: string,
}

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: IFunctionArgs
    },
    scripts?: {
        [_: string]: IFunctionArgs
    }
    queries?: {
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
    queryFunctionName?: string,
    isDataSourceQuery?: boolean,
    layoutViewState?: string,
}

export const CampaignFieldTypes = {
  String: DG.TYPE.STRING,
  Number: DG.TYPE.FLOAT,
  Boolean: DG.TYPE.BOOL,
  Date: DG.TYPE.DATE_TIME,
  Molecule: DG.SEMTYPE.MOLECULE,
} as const;

export type HitTriageCampaignFieldType = keyof typeof CampaignFieldTypes;

export type HitTriageCampaignField =
    {name: string, type: HitTriageCampaignFieldType, semtype?: string, required?: boolean};

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
    scripts?: HitTriageTemplateScript[],
    queries?: HitTriageTemplateScript[],
};

export type HitTriageTemplateSubmit = {
    fName: string,
    package: string
};

export type HitTriageCampaignStatus = string;

export type HitTriageCampaign = {
    name: string,
    templateName: string,
    template?: HitDesignTemplate,
    savePath?: string,
    status: HitTriageCampaignStatus,
    createDate: string,
    campaignFields: {[key: string]: any},
    filters: {[key: string]: any}[],
    ingest: HitTriageTemplateIngest,
    columnSemTypes?: {[key: string]: string},
    rowCount?: number,
    filteredRowCount?: number,
    layout?: string,
    columnTypes?: {[key: string]: string},
    version?: number,
    permissions?: TriagePermissions,
    authorUserId?: string
};

export type TriagePermissions = {
    // will store the ids of UserGroups that have access to the campaign
    edit: string[],
    view: string[],
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

export type HitDesignTemplate = Omit<HitTriageTemplate, 'dataSourceType' | 'queryFunctionName'> &
    {stages: string[]};

export type HitDesignCampaign = Omit<HitTriageCampaign, 'filters' | 'ingest'> & {tilesViewerFormSketch?: string};

// todo: probably add some more stuff
export type PeptiHitTemplate = HitDesignTemplate & {toAtomiLevelProps?: {[key: string]: any}}

export type ComputeFunctions = {
    functions: DG.Func[],
    scripts: DG.Script[],
    queries: DG.DataQuery[],
};
