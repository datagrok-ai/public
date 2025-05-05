import {Subsector} from '../sparklines/piechart';

export enum SectorType {
  SECTOR = 'sector',
  SUBSECTOR = 'subsector'
}

export interface VlaaivisColumnMetadata extends Subsector {
  groupName?: string;
  sectorColor?: string;
}

export const TAGS = {
  VLAAIVIS_METADATA: '.vlaaivis-metadata',
};

export const CONSTANTS = {
  SECTOR_COLOR_PROPERTY: 'sectorColor',
  LOWER_BOUND: 'lowerBound',
  UPPER_BOUND: 'upperBound',
  VLAAIVIS: 'VlaaiVis'
};

export const DEFAULTS = {
  WEIGHT: '0.2'
};

export const groupProps = [
  {
    'property': {
      'name': 'sectorColor',
      'inputType': 'Color',
      'description': 'Defines the color used for sector'
    },
    'object': {
      'sectorColor': '#e66465'
    }
  }
];

export const generalProps = [
  {
    'property': {
      'name': 'lowerBound',
      'inputType': 'Float',
      'description': 'Defines the minimum acceptable value within the TPP'
    },
    'object': {
      'lowerBound': 0.8
    }
  },
  {
    'property': {
      'name': 'upperBound',
      'inputType': 'Float',
      'description': 'Defines the maximum acceptable value within the TPP'
    },
    'object': {
      'upperBound': 0.9
    }
  }
];

export const subGroupProps = [
  {
    'property': {
      'name': 'min',
      'inputType': 'Float',
      'enabled': false,
      'description': 'min value of the property (optional; used for editing the line)'
    },
    'object': {}
  },
  {
    'property': {
      'name': 'max',
      'inputType': 'Float',
      'enabled': false,
      'description': 'max value of the property (optional; used for editing the line)'
    },
    'object': {}
  },
  {
    'property': {
      'name': 'weight',
      'inputType': 'Slider',
      'min': 0,
      'max': 1,
      'description': 'Determines the thickness of the wedge'
    },
    'object': {
      'weight': 0.2
    }
  }
];

export const defaultGeneralProps = generalProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record <string, any>);

export const defaultGroupProps = groupProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record <string, any>);
