export interface Sample {
  [key: string]: number | string | undefined;
}

export interface DescriptorStats {
  name: string;
  meanGood: number;
  meanBad: number;
  stdGood: number;
  stdBad: number;
  weight: number;
  t: number;
  p: number;
}
