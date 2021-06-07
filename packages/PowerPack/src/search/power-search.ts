/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from "../widgets/web-widget";

// Power Search: community-curated, template-based, widget-driven search engine

interface Template {
  template: string;
  url: string;
  regexp?: RegExp;   // cached regexp for the template
}

interface Card {
  id: string;
  name: string;
  templates: Template[];
}

export function powerSearch(s: string, host: HTMLDivElement): void {
  for (let p of templates)
    for (let t of p.templates) {
      t.regexp ??= new RegExp(t.template);
      let matches = t.regexp.exec(s);
      let url = t.url;

      if (matches !== null) {
        for (let i = 1; i < matches.length; i++)
          url = url.replace('${' + i + '}', matches[i]);

        let widget = new WebWidget({
            src: url,
            width: '100%',
            height: '500px'
          });

        host.appendChild(widget.root);
      }
    }
}

const semTypes = [
  {
    name: 'CHEMBL_ID',
    description: 'ChEMBL compound identifier',
    template: '(CHEMBL[0-9]+)'
  }
];

const templates: Card[] = [
  {
    id: 'chembl-id-report-card',
    name: 'Report card',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/name_and_classification'
    }]
  },
  {
    id: 'chembl-id-representations',
    name: 'Representations',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/representations'
    }]
  },
  {
    id: 'chembl-id-alternative-forms',
    name: 'Alternative forms',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/alternate_forms'
    }]
  },
  {
    id: 'chembl-id-calculated-properties',
    name: 'Calculated properties',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/calculated_properties'
    }]
  },  {
    id: 'chembl-id-cross-refs',
    name: 'Cross references',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/cross_refs'
    }]
  },
]

const widgetTemplates = [

];

// <object data="https://www.ebi.ac.uk/chembl/embed/#compound_report_card/CHEMBL1193654/name_and_classification" width="100%" height="100%"></object>