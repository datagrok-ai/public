import { Plate } from "./plate";

export interface IPlateWellValidator {
  name: string;
  description: string;
  validate(plate: Plate, row: number, col: number): string | null;
}

export const plateWellValidators: IPlateWellValidator[] = [
  {
    name: 'Concentration and volume go together',
    description: 'Concentration and volume should either be both present or both absent',
    validate: (plate: Plate, row: number, col: number) => {
      if (!plate.data.col('concentration') || !plate.data.col('volume')) 
        return null;

      const concentration = plate.get('concentration', row, col);
      const volume = plate.get('volume', row, col);

      if ((concentration && !volume) || (!concentration && volume))
        return 'Concentration and volume should either be both present or both absent';

      return null;
    },
  },

  {
    name: 'Concentration is not negative',
    description: 'Concentration should not be negative',
    validate: (plate: Plate, row: number, col: number) => {
      if (!plate.data.col('concentration'))
        return null;

      const concentration = plate.get('concentration', row, col);
      if (concentration != null && concentration < 0)
        return 'Concentration should not be negative';

      return null;  
    },
  },
];