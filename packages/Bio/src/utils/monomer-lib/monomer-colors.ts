import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Helm} from '../helm-to-molfile/converter/helm';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';

/**
 * MonomerColors class from HelmWebEditor for natural monomers
 */
export const naturalMonomerColors = {
  [HelmTypes.BASE]: {
    A: "#A0A0FF",
    G: "#FF7070",
    T: "#A0FFA0",
    C: "#FF8C4B",
    U: "#FF8080"
  },

  [HelmTypes.NUCLEOTIDE]: {
    A: "#A0A0FF",
    G: "#FF7070",
    T: "#A0FFA0",
    C: "#FF8C4B",
    U: "#FF8080"
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
    A: "#C8C8C8",
    R: "#145AFF",
    N: "#00DCDC",
    D: "#E60A0A",
    C: "#E6E600",
    E: "#00DCDC",
    Q: "#E60A0A",
    G: "#EBEBEB",
    H: "#8282D2",
    I: "#0F820F",
    L: "#0F820F",
    K: "#145AFF",
    M: "#E6E600",
    F: "#3232AA",
    P: "#DC9682",
    S: "#FA9600",
    T: "#FA9600",
    W: "#B45AB4",
    Y: "#3232AA",
    V: "#0F820F"
  },

  [HelmTypes.CHEM]: {
    R: "#eeeeee",
  },

  [HelmTypes.BLOB]: {
    B: "#999999",
    G: "#e2e2e2"
  }
};

