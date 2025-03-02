import {Plate} from "./plate";
import {BmgPherastar96PlateReader} from "./readers/bmg-pherastar-96";
import {DelfiaEnvision96MaxRluPlateReader} from "./readers/delfia-envision-96-max-rlu";
import {SpectramaxPlateReader} from "./readers/spectramax";
import {DelfiaEnvisionByRowPlateReader} from "./readers/delfia-envision-by-row";

export interface IPlateReader {
  isApplicable(s: string): boolean;
  read(s: string): Plate;
}


export abstract class PlateReader {
  static readers: IPlateReader[] = [
    new DelfiaEnvision96MaxRluPlateReader(),
    new DelfiaEnvisionByRowPlateReader(),
    new BmgPherastar96PlateReader(),
    new SpectramaxPlateReader(),
  ];

  static getReader(s: string): IPlateReader | null {
    return this.readers.find((r) => r.isApplicable(s)) ?? null;
  }

  static read(s: string): Plate | null {
    return this.getReader(s)?.read(s) ?? null;
  }
}
