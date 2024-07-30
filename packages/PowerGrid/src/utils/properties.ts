export const groupProps = [
  {
    "property": {
      "name": "sectorColor",
      "inputType": "Color",
      "description": "Defines the color used for sector"
    },
    "object": {
      "sectorColor": "#e66465"
    }
  }
];

export const generalProps = [
  {
    "property": {
      "name": "lowerBound",
      "inputType": "Float",
      "description": "Defines the minimum acceptable value within the TPP"
    },
    "object": {
      "lowerBound": 0.8
    }
  },
  {
    "property": {
      "name": "upperBound",
      "inputType": "Float",
      "description": "Defines the maximum acceptable value within the TPP"
    },
    "object": {
      "upperBound": 0.9
    }
  }
];

export const subGroupProps = [
  {
    "property": {
      "name": "low",
      "inputType": "Float",
      "description": "Sets the minimum value for a sub-sector"
    },
    "object": {
      "low": 0
    }
  },
  {
    "property": {
      "name": "high",
      "inputType": "Float",
      "description": "Sets the maximum value for a sub-sector"
    },
    "object": {
      "high": 1
    }
  },
  {
    "property": {
      "name": "weight",
      "inputType": "Float",
      "description": "Determines the thickness of the wedge"
    },
    "object": {
      "weight": 0.2
    }
  }
];