import {Plate} from "./plate";
import {BmgPherastar96PlateReader} from "./readers/bmg-pherastar-96";
import {DelfiaEnvision96MaxRluPlateReader} from "./readers/delfia-envision-96-max-rlu";

export interface IPlateReader {
  isApplicable(s: string): boolean;
  read(s: string): Plate;
}


export abstract class PlateReader {
  static readers: IPlateReader[] = [
    new DelfiaEnvision96MaxRluPlateReader(),
    new BmgPherastar96PlateReader(),
  ];

  static getReader(s: string): IPlateReader | null {
    return this.readers.find((r) => r.isApplicable(s)) ?? null;
  }

  static read(s: string): Plate | null {
    return this.getReader(s)?.read(s) ?? null;
  }
}
