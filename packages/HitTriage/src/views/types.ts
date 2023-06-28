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

export type IChemProperty = {
   propertyName: string
} & Descriptor;

export type IChemPropertyGroup = {
    name: ChemPropNames,
    values: IChemProperty[],
}

export type IChemPropertyGroupMap = {
    toxRisks: IChemPropertyGroup,
    structuralAlerts: IChemPropertyGroup,
    chemProperties: IChemPropertyGroup,
}

export type ChemPropNames = 'Descriptors' | 'Toxicity Risks' | 'Structural Alerts' | 'Chemical Properties';

export type ChemFunctionType = (t: DG.DataFrame, col: string, props: string[]) => Promise<any>;
