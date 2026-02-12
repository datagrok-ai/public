import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {studies} from '../utils/app-utils';
import {SUBJECT_ID, DOMAIN, PLANNED_TRT_ARM, ACT_TRT_ARM} from '../constants/columns-constants';
import {StudyTableViewParams} from '../utils/views-creation-utils';
import {createAllMeasurementsDf} from '../data-preparation/data-preparation';
import {createSubjectProfileView} from './subject-profile-view';
import {awaitCheck} from '@datagrok-libraries/test/src/test';
import {subjectClicked$, PADDING_X, MARKER_RADIUS, getArmColorStr} from '../utils/combined-measurements-cell-renderer';
import * as PinnedUtils from '@datagrok-libraries/gridext/src/pinned/PinnedUtils';

type ArmEntry = {values: number[], subjects: string[]};
type ArmValues = {[arm: string]: ArmEntry};
type DayMap = {[day: number]: ArmValues};
type GroupData = {domain: string, category: string, test: string, units: string, days: DayMap};

export function createEventsView(studyId: string): StudyTableViewParams {
  const allMeasurements = createAllMeasurementsDf(studyId);
  if (!allMeasurements)
    return {df: DG.DataFrame.create()};

  const dm = studies[studyId].domains.dm;
  if (dm) {
    const columnsFromDm = dm.columns.names().filter((it) => it !== SUBJECT_ID);
    grok.data.joinTables(allMeasurements, dm, [SUBJECT_ID], [SUBJECT_ID], null,
      columnsFromDm, DG.JOIN_TYPE.LEFT, true);
  }

  const intVisitDay = (allMeasurements as DG.DataFrame).col('visit_day')!.convertTo( DG.TYPE.INT);
  (allMeasurements as DG.DataFrame).columns.remove('visit_day');
  (allMeasurements as DG.DataFrame).columns.add(intVisitDay);

  const armColName = allMeasurements.col(PLANNED_TRT_ARM)
    ? PLANNED_TRT_ARM
    : (allMeasurements.col(ACT_TRT_ARM) ? ACT_TRT_ARM : null);

  const domainCol = allMeasurements.col(DOMAIN)!;
  const catCol = allMeasurements.col('category')!;
  const testCol = allMeasurements.col('test')!;
  const unitsCol = allMeasurements.col('units')!;
  const resultCol = allMeasurements.col('result')!;
  const visitDayCol = allMeasurements.col('visit_day')!;
  const armCol = armColName ? allMeasurements.col(armColName) : null;
  const subjectCol = allMeasurements.col(SUBJECT_ID)!;

  const uniqueDays = new Set<number>();
  const uniqueArms = new Set<string>();
  const uniqueSubjects = new Set<string>();
  const groupMap = new Map<string, GroupData>();

  for (let i = 0; i < allMeasurements.rowCount; i++) {
    if (resultCol.isNone(i) || visitDayCol.isNone(i))
      continue;

    const domain = domainCol.isNone(i) ? '' : domainCol.get(i);
    const category = catCol.isNone(i) ? '' : catCol.get(i);
    const test = testCol.isNone(i) ? '' : testCol.get(i);
    const units = unitsCol.isNone(i) ? '' : unitsCol.get(i);
    const day = visitDayCol.get(i);
    const result = resultCol.get(i);
    const arm = armCol ? (armCol.isNone(i) ? 'Unknown' : armCol.get(i)) : 'Unknown';
    const subject = subjectCol.isNone(i) ? '' : subjectCol.get(i);
    const groupKey = `${domain}|${category}|${test}|${units}`;

    uniqueDays.add(day);
    uniqueArms.add(arm);
    if (subject)
      uniqueSubjects.add(subject);

    if (!groupMap.has(groupKey))
      groupMap.set(groupKey, {domain, category, test, units, days: {}});

    const group = groupMap.get(groupKey)!;
    if (!group.days[day])
      group.days[day] = {};
    if (!group.days[day][arm])
      group.days[day][arm] = {values: [], subjects: []};
    group.days[day][arm].values.push(result);
    group.days[day][arm].subjects.push(subject);
  }

  const sortedDays = [...uniqueDays].sort((a, b) => a - b);
  const sortedArms = [...uniqueArms].sort();

  const rowCount = groupMap.size;
  const pivotedDf = DG.DataFrame.create(rowCount);
  pivotedDf.name = 'Events';

  const pDomainCol = pivotedDf.columns.addNewString(DOMAIN);
  const pCatCol = pivotedDf.columns.addNewString('category');
  const pTestCol = pivotedDf.columns.addNewString('test');
  const pUnitsCol = pivotedDf.columns.addNewString('units');

  const yAxisCol = pivotedDf.columns.addNewString('Y axis');
  yAxisCol.semType = 'y-axis';

  const dayColumns: DG.Column[] = [];
  for (const day of sortedDays) {
    const col = pivotedDf.columns.addNewString(`Day ${day}`);
    col.semType = 'combined-measurements';
    dayColumns.push(col);
  }

  const rowMinCol = pivotedDf.columns.addNewFloat('~rowMin');
  const rowMaxCol = pivotedDf.columns.addNewFloat('~rowMax');

  let rowIdx = 0;
  for (const [_key, group] of groupMap) {
    pDomainCol.set(rowIdx, group.domain);
    pCatCol.set(rowIdx, group.category);
    pTestCol.set(rowIdx, group.test);
    pUnitsCol.set(rowIdx, group.units);

    let rowMin = Infinity;
    let rowMax = -Infinity;
    for (let d = 0; d < sortedDays.length; d++) {
      const dayData = group.days[sortedDays[d]];
      if (dayData) {
        dayColumns[d].set(rowIdx, JSON.stringify(dayData));
        for (const arm of Object.keys(dayData)) {
          for (const v of dayData[arm].values) {
            if (v < rowMin) rowMin = v;
            if (v > rowMax) rowMax = v;
          }
        }
      }
    }
    const finalMin = rowMin === Infinity ? 0 : rowMin;
    const finalMax = rowMax === -Infinity ? 0 : rowMax;
    rowMinCol.set(rowIdx, finalMin);
    rowMaxCol.set(rowIdx, finalMax);
    yAxisCol.set(rowIdx, JSON.stringify({min: finalMin, max: finalMax}));
    rowIdx++;
  }

  const sortedSubjects = [...uniqueSubjects].sort();
  pivotedDf.setTag('armOrder', JSON.stringify(sortedArms));
  pivotedDf.setTag('subjectOrder', JSON.stringify(sortedSubjects));

  const onTableViewAdded = async (tableView: DG.TableView) => {
    await awaitCheck(() => tableView.grid !== null, '', 1000);

    const grid = tableView.grid;
    const metaCols = [DOMAIN, 'category', 'test', 'units'];
    const dayColNames = sortedDays.map((d) => `Day ${d}`);
    grid.columns.setOrder(metaCols.concat(['Y axis']).concat(dayColNames));

    const yAxisGridCol = grid.col('Y axis');
    if (yAxisGridCol)
      yAxisGridCol.width = 50;

    for (const dayColName of dayColNames) {
      const gridCol = grid.col(dayColName);
      if (gridCol)
        gridCol.width = 120;
    }

    const minGridCol = grid.col('~rowMin');
    const maxGridCol = grid.col('~rowMax');
    if (minGridCol) minGridCol.visible = false;
    if (maxGridCol) maxGridCol.visible = false;

    grid.setOptions({rowHeight: 60, colHeaderHeight: 50});

    const dayColNameSet = new Set(dayColNames);

    const HEADER_HIT_RADIUS = 6;

    grid.onCellRendered.subscribe((args) => {
      const cell = args.cell;
      if (!cell.isColHeader)
        return;
      const colName = cell.tableColumn?.name;
      if (!colName || !dayColNameSet.has(colName))
        return;

      const g = args.g;
      const bounds = args.bounds;
      const markerY = bounds.y + bounds.height - MARKER_RADIUS - 4;
      const armWidth = (bounds.width - 2 * PADDING_X) / sortedArms.length;

      g.lineWidth = 1.5;
      for (let i = 0; i < sortedArms.length; i++) {
        const cx = bounds.x + PADDING_X + i * armWidth + armWidth / 2;
        g.strokeStyle = getArmColorStr(i);
        g.beginPath();
        g.arc(cx, markerY, MARKER_RADIUS, 0, 2 * Math.PI);
        g.stroke();
      }
    });

    grid.root.addEventListener('mousemove', (e: MouseEvent) => {
      const cell = grid.hitTest(e.offsetX, e.offsetY);
      if (!cell || !cell.isColHeader) return;

      const colName = cell.tableColumn?.name;
      if (!colName || !dayColNameSet.has(colName)) return;

      const bounds = cell.bounds;
      const markerY = bounds.y + bounds.height - MARKER_RADIUS - 4;
      const armWidth = (bounds.width - 2 * PADDING_X) / sortedArms.length;

      for (let i = 0; i < sortedArms.length; i++) {
        const cx = bounds.x + PADDING_X + i * armWidth + armWidth / 2;
        const dx = e.offsetX - cx;
        const dy = e.offsetY - markerY;
        if (Math.sqrt(dx * dx + dy * dy) < HEADER_HIT_RADIUS) {
          const el = ui.divText(sortedArms[i]);
          el.style.color = getArmColorStr(i);
          ui.tooltip.show(el, e.clientX + 16, e.clientY + 16);
          return;
        }
      }
    });

    for (const colName of metaCols.concat(['Y axis'])) {
      const gridCol = grid.col(colName);
      if (gridCol)
        PinnedUtils.addPinnedColumn(gridCol);
    }

    const subjectInput = ui.input.choice('Subject', {
      value: '',
      items: ['', ...sortedSubjects],
      onValueChanged: (value) => {
        pivotedDf.setTag('selectedSubject', value);
        grid.invalidate();
        if (value) {
          const profileView = createSubjectProfileView(studyId, allMeasurements, value);
          grok.shell.o = profileView;
        }
      },
    });
    const tooltipContent = ui.divV([]);
    const legendDiv = ui.div([], {classes: 'preclinical-case-legend'});
    for (let i = 0; i < sortedArms.length; i++) {
      const color = DG.Color.getCategoricalColor(i);
      const r = (color >> 16) & 0xFF;
      const g = (color >> 8) & 0xFF;
      const b = color & 0xFF;
      const colorStr = `rgb(${r},${g},${b})`;

      const dot = ui.div([], {classes: 'preclinical-case-legend-dot', style: {borderColor: colorStr}});
      legendDiv.appendChild(dot);

      const tooltipDot = ui.div([], {classes: 'preclinical-case-legend-tooltip-dot', style: {borderColor: colorStr}});
      const tooltipLabel = ui.div([sortedArms[i]], {classes: 'preclinical-case-legend-tooltip-label', style: {color: colorStr}});
      const tooltipRow = ui.divH([tooltipDot, tooltipLabel], {classes: 'preclinical-case-legend-tooltip-row'});
      tooltipContent.appendChild(tooltipRow);
    }
    ui.tooltip.bind(legendDiv, () => tooltipContent);

    tableView.setRibbonPanels([[subjectInput.root, legendDiv]]);

    subjectClicked$.subscribe((subject) => {
      subjectInput.value = subject;
      grid.invalidate();
    });

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL,
      column: DOMAIN,
      columnName: DOMAIN,
    });
  };

  return {df: pivotedDf, onTableViewAddedFunc: onTableViewAdded};
}