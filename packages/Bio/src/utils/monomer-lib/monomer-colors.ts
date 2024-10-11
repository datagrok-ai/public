import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Helm} from '../helm-to-molfile/converter/helm';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';

/**
 * MonomerColors class from HelmWebEditor for natural monomers
 */
export const naturalMonomerColors = {
  [HelmTypes.BASE]: {
    // Chromatogram palette  // HELMWebEditor monomerColors
    A: "#20E040",            // "#A0A0FF",
    G: "#040404",            // "#FF7070",
    T: "#FF8080",            // "#A0FFA0",
    C: "#2060FF",            // "#FF8C4B",
    U: "#FF8080",            // "#FF8080"
  },

  [HelmTypes.NUCLEOTIDE]: {
    // Chromatogram palette  // HELMWebEditor monomerColors
    A: "#20E040",            // "#A0A0FF",
    G: "#040404",            // "#FF7070",
    T: "#FF8080",            // "#A0FFA0",
    C: "#2060FF",            // "#FF8C4B",
    U: "#FF8080",            // "#FF8080"
  },

  [HelmTypes.LINKER]: {
    P: "#9aa5e1",
    p: "#9aa5e1"
  },

  [HelmTypes.SUGAR]: {
    R: "#7a85c1",
    r: "#7a85c1",
    // TODO: deoxyribose
  },

  [HelmTypes.AA]: {
    // GrokGroups palette    // HELMWebEditor monomerColors
    A: "rgb(44,160,44)",     // "#C8C8C8",
    R: "rgb(23,190,207)",    // "#145AFF",
    N: "rgb(235,137,70)",    // "#00DCDC",
    D: "rgb(31,119,180)",    // "#E60A0A",
    C: "rgb(188,189,34)",    // "#E6E600",
    E: "rgb(31, 120, 150)",  // "#00DCDC",
    Q: "rgb(205, 111, 71)",  // "#E60A0A",
    G: "rgb(214,39,40)",     // "#EBEBEB",
    H: "rgb(158,218,229)",   // "#8282D2",
    I: "rgb(23,103,57)",     // "#0F820F",
    L: "rgb(30,110,96)",     // "#0F820F",
    K: "rgb(108, 218, 229)", //"#145AFF",
    M: "rgb(60,131,95)",     // "#E6E600",
    F: "rgb(24,110,79)",     // "#3232AA",
    P: "rgb(255,152,150)",   // "#DC9682",
    S: "rgb(255,187,120)",   // "#FA9600",
    T: "rgb(245,167,100)",   // "#FA9600",
    W: "rgb(182, 223, 138)", // "#B45AB4",
    Y: "rgb(152,223,138)",   // "#3232AA",
    V: "rgb(74,160,74)",     // "#0F820F",
  },

  [HelmTypes.CHEM]: {
    R: "#eeeeee",
  },

  [HelmTypes.BLOB]: {
    B: "#999999",
    G: "#e2e2e2"
  }
};

