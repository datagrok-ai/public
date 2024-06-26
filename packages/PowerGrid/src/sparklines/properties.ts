export const groupProps = [
  {
    "property": {
      "name": "sectorColor",
      "inputType": "Color"
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
      "inputType": "Float"
    },
    "object": {
      "lowerBound": 0.8
    }
  },
  {
    "property": {
      "name": "upperBound",
      "inputType": "Float"
    },
    "object": {
      "upperBound": 0.9
    }
  }
];

export const subGroupProps = [
  {
    "property": {
      "name": "lowThreshold",
      "inputType": "Float"
    },
    "object": {
      "lowThreshold": 0 //use min column value
    }
  },
  {
    "property": {
      "name": "highThreshold",
      "inputType": "Float"
    },
    "object": {
      "highThreshold": 1 //use max column value
    }
  },
  {
    "property": {
      "name": "weight",
      "inputType": "Float"
    },
    "object": {
      "weight": 0.2 //generate here numbers too?
    }
  }
]