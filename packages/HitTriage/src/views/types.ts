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
