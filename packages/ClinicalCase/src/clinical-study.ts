/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as meta from './sdtm-meta';

export class ClinicalDomains {
  ae: DG.DataFrame = null;
  ce: DG.DataFrame = null;
  ds: DG.DataFrame = null;
  dv: DG.DataFrame = null;
  ho: DG.DataFrame = null;
  mh: DG.DataFrame = null;
  cv: DG.DataFrame = null;
  da: DG.DataFrame = null;
  dd: DG.DataFrame = null;
  eg: DG.DataFrame = null;
  fa: DG.DataFrame = null;
  ft: DG.DataFrame = null;
  ie: DG.DataFrame = null;
  is: DG.DataFrame = null;
  lb: DG.DataFrame = null;
  mb: DG.DataFrame = null;
  mi: DG.DataFrame = null;
  mk: DG.DataFrame = null;
  mo: DG.DataFrame = null;
  ms: DG.DataFrame = null;
  nv: DG.DataFrame = null;
  oe: DG.DataFrame = null;
  pc: DG.DataFrame = null;
  pe: DG.DataFrame = null;
  pp: DG.DataFrame = null;
  qs: DG.DataFrame = null;
  re: DG.DataFrame = null;
  rp: DG.DataFrame = null;
  rs: DG.DataFrame = null;
  sc: DG.DataFrame = null;
  sr: DG.DataFrame = null;
  ss: DG.DataFrame = null;
  tr: DG.DataFrame = null;
  tu: DG.DataFrame = null;
  ur: DG.DataFrame = null;
  vs: DG.DataFrame = null;
  ag: DG.DataFrame = null;
  cm: DG.DataFrame = null;
  ec: DG.DataFrame = null;
  ex: DG.DataFrame = null;
  ml: DG.DataFrame = null;
  pr: DG.DataFrame = null;
  su: DG.DataFrame = null;
  co: DG.DataFrame = null;
  dm: DG.DataFrame = null;
  se: DG.DataFrame = null;
  sm: DG.DataFrame = null;
  sv: DG.DataFrame = null;
  ta: DG.DataFrame = null;
  td: DG.DataFrame = null;
  te: DG.DataFrame = null;
  ti: DG.DataFrame = null;
  tm: DG.DataFrame = null;
  ts: DG.DataFrame = null;
  tv: DG.DataFrame = null;
}

export class ClinicalStudy {
  name: string;
  description: string;
  domains: ClinicalDomains = new ClinicalDomains();
  subjectsCount: number;
  sitesCount: number;

  initFromWorkspace(): void {
    for (let t of grok.shell.tables) {
      if (t.name.toLowerCase() in this.domains)
        this.domains[t.name.toLowerCase()] = t;
    }

    this.subjectsCount = this.domains.dm.rowCount;
    this.sitesCount = this.domains.dm.col('siteid').stats.uniqueCount;
    this.name = this.domains.dm.col('studyid').get(0);
  }
}

export let study: ClinicalStudy = new ClinicalStudy();
