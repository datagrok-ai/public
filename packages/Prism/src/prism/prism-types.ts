export interface PrismFile {
  sheets: PrismSheet[];
  analyses: PrismAnalysis[];
  sheetNames: Map<string, string>;
}

export interface PrismSheet {
  id: string;
  title: string;
  dataFormat: string;
  replicatesCount: number;
  xTitle: string | null;
  rowTitlesPresent: boolean;
  yColumnTitles: string[];
  data: string[][];
}

export interface PrismAnalysis {
  id: string;
  title: string;
  resultData: string[][];
  columnTitles: string[];
}

export interface DocumentJson {
  sheets: {
    data?: string[];
    analyses?: string[];
    graphs?: string[];
    info?: string[];
    layouts?: string[];
  };
  sheetAttributesMap: Record<string, any>;
}

export interface SheetJson {
  title: string;
  table: {
    uid: string;
    dataFormat: string;
    replicatesCount?: number;
    dataSets: string[];
    xDataSet?: string;
    rowTitlesDataSet?: string;
  };
  preferredResultSheet?: string;
}

export interface DatasetJson {
  uid: string;
  title?: string;
  dataSheet?: string;
}
