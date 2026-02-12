import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Dayjs} from 'dayjs';
import {SUMMARY_VIEW_NAME} from '../constants/view-names-constants';
import {StudyConfig} from '../types/types';
import {openStudy, studies} from './app-utils';

export function createInitialSatistics(clinicalCaseNode: DG.TreeViewGroup, configs: StudyConfig[]): HTMLDivElement {
  const tableDiv = ui.div('', {style: {position: 'relative', paddingTop: '15px'}});

  const statsObj: {[key: string]: string}[] = [];
  for (const config of configs) {
    const item: {[key: string]: any} = {};
    for (const prop of Object.keys(config)) {
      if (prop !== 'other') {
        if (prop === 'startDate' || prop === 'endDate')
          item[prop] = (config[prop] as Dayjs).toString();
        else
          item[prop] = (config as any)[prop] ?? '';
      }
    }
    statsObj.push(item);
  }

  const table = ui.table(statsObj, (item) => {
    return [
      ui.link(item.name, async () => {
        openStudy(clinicalCaseNode, item.name, SUMMARY_VIEW_NAME);
      }, 'Open study'),
      item.protocol,
      item.totalSubjects,
      item.startDate,
      item.endDate,
    ];
  },
  ['Study', 'Protocol', 'Total subjects', 'Start date', 'End date']);

  tableDiv.append(table);
  return tableDiv;
}
