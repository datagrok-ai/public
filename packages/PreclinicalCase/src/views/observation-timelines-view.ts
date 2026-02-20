import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {studies} from '../utils/app-utils';
import {SUBJECT_ID, DOMAIN, PLANNED_TRT_ARM, ACT_TRT_ARM} from '../constants/columns-constants';
import {StudyTableViewParams} from '../utils/views-creation-utils';
import {restoreBrowsePanelOnRemoval} from '../utils/utils';
import {createAllMeasurementsDf} from '../data-preparation/data-preparation';
import {createSubjectProfileView} from './subject-profile-view';
import {awaitCheck} from '@datagrok-libraries/test/src/test';
import {subjectClicked$, PADDING_X, PADDING_Y, getArmColorStr, getEnabledArms} from '../utils/combined-measurements-cell-renderer';
import * as PinnedUtils from '@datagrok-libraries/gridext/src/pinned/PinnedUtils';
import {Subscription} from 'rxjs';

type ArmEntry = {values: number[], subjects: string[]};
type ArmValues = {[arm: string]: ArmEntry};
type DayMap = {[day: number]: ArmValues};
type GroupData = {domain: string, category: string, test: string, units: string, days: DayMap};

let subjectClickedSub: Subscription | null = null;

export function createObservationTimelinesView(studyId: string): StudyTableViewParams {
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
  pivotedDf.name = 'Observation timelines';

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

  const dayColDataMasks: DG.BitSet[] = [];
  for (const col of dayColumns) {
    const mask = DG.BitSet.create(rowCount);
    for (let r = 0; r < rowCount; r++) {
      if (!col.isNone(r))
        mask.set(r, true);
    }
    dayColDataMasks.push(mask);
  }

  const armCodeToName = studies[studyId].armCodeToName;
  const getArmCode = (arm: string) => Object.keys(armCodeToName).find((c) => armCodeToName[c] === arm) ?? arm;

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

    grid.setOptions({rowHeight: 60, colHeaderHeight: 100});

    const dayColNameSet = new Set(dayColNames);

    grid.onCellRendered.subscribe((args) => {
      const cell = args.cell;
      if (!cell.isColHeader)
        return;
      const colName = cell.tableColumn?.name;
      if (!colName || !dayColNameSet.has(colName))
        return;

      const g = args.g;
      const bounds = args.bounds;
      const curEnabledArms = getEnabledArms(pivotedDf);
      if (curEnabledArms.length === 0)
        return;
      const armWidth = (bounds.width - 2 * PADDING_X) / curEnabledArms.length;

      const maxLabelHeight = 30;
      g.font = '10px sans-serif';
      g.textBaseline = 'middle';
      for (let i = 0; i < curEnabledArms.length; i++) {
        const globalIdx = sortedArms.indexOf(curEnabledArms[i]);
        const cx = bounds.x + PADDING_X + i * armWidth + armWidth / 2;
        let armCode = getArmCode(curEnabledArms[i]);
        while (armCode.length > 1 && g.measureText(armCode).width > maxLabelHeight)
          armCode = armCode.slice(0, -1);
        if (armCode !== getArmCode(curEnabledArms[i]))
          armCode += '\u2026';
        g.save();
        g.translate(cx, bounds.y + bounds.height - 4);
        g.rotate(-Math.PI / 2);
        g.fillStyle = getArmColorStr(globalIdx >= 0 ? globalIdx : i);
        g.fillText(armCode, 0, 0);
        g.restore();
      }
    });

    grid.root.addEventListener('mousemove', (e: MouseEvent) => {
      if (e.target !== grid.canvas && e.target !== grid.overlay)
        return;
      const cell = grid.hitTest(e.offsetX, e.offsetY);
      if (!cell || !cell.isColHeader) return;

      const colName = cell.tableColumn?.name;
      if (!colName || !dayColNameSet.has(colName)) return;

      const bounds = cell.bounds;
      const curEnabledArms = getEnabledArms(pivotedDf);
      if (curEnabledArms.length === 0) return;
      const armWidth = (bounds.width - 2 * PADDING_X) / curEnabledArms.length;

      for (let i = 0; i < curEnabledArms.length; i++) {
        const bandStart = bounds.x + PADDING_X + i * armWidth;
        if (e.offsetX >= bandStart && e.offsetX < bandStart + armWidth) {
          const globalIdx = sortedArms.indexOf(curEnabledArms[i]);
          const el = ui.divText(curEnabledArms[i]);
          el.style.color = getArmColorStr(globalIdx >= 0 ? globalIdx : i);
          ui.tooltip.show(el, e.clientX + 16, e.clientY + 16);
          return;
        }
      }
    });

    for (const colName of metaCols) {
      const gridCol = grid.col(colName);
      if (!gridCol) continue;
      if (colName === 'test')
        gridCol.width = 80;
      else
        gridCol.visible = false;
    }

    for (const colName of ['test', 'Y axis']) {
      const pinCol = grid.col(colName);
      if (pinCol)
        PinnedUtils.addPinnedColumn(pinCol);
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
    const disabledArms = new Set<string>();
    const baseArmWidth = (120 - 2 * PADDING_X) / sortedArms.length;
    const legendDiv = ui.div([], {classes: 'preclinical-case-legend'});
    for (let i = 0; i < sortedArms.length; i++) {
      const color = DG.Color.getCategoricalColor(i);
      const r = (color >> 16) & 0xFF;
      const g = (color >> 8) & 0xFF;
      const b = color & 0xFF;
      const colorStr = `rgb(${r},${g},${b})`;

      const armCode = getArmCode(sortedArms[i]);
      const label = ui.div([armCode], {classes: 'preclinical-case-legend-label', style: {color: colorStr}});
      label.style.cursor = 'pointer';
      ui.tooltip.bind(label, () => {
        const tip = ui.divText(sortedArms[i]);
        tip.style.color = disabledArms.has(sortedArms[i]) ? '#aaa' : colorStr;
        return tip;
      });
      label.addEventListener('click', () => {
        if (disabledArms.has(sortedArms[i])) {
          disabledArms.delete(sortedArms[i]);
          label.style.color = colorStr;
          label.style.opacity = '';
        } else {
          disabledArms.add(sortedArms[i]);
          label.style.color = '#ccc';
          label.style.opacity = '0.5';
        }
        pivotedDf.setTag('disabledArms', JSON.stringify([...disabledArms]));
        const enabledCount = sortedArms.length - disabledArms.size;
        const newWidth = Math.round(2 * PADDING_X + enabledCount * baseArmWidth);
        for (const dayColName of dayColNames) {
          const gc = grid.col(dayColName);
          if (gc && gc.visible)
            gc.width = newWidth;
        }
        grid.invalidate();
      });
      legendDiv.appendChild(label);
    }

    const connectPointsInput = ui.input.bool('Connect subject points', {value: false, onValueChanged: () => grid.invalidate()});
    tableView.setRibbonPanels([[subjectInput.root, connectPointsInput.root, legendDiv]]);

    subjectClickedSub?.unsubscribe();
    subjectClickedSub = subjectClicked$.subscribe((subject) => {
      subjectInput.value = subject;
      grid.invalidate();
    });

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    for (const colName of metaCols) {
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: colName,
        columnName: colName,
      });
    }

    pivotedDf.onFilterChanged.subscribe(() => {
      const filter = pivotedDf.filter;
      const visible = new Uint8Array(dayColumns.length);
      for (let r = filter.findNext(-1, true); r !== -1; r = filter.findNext(r, true)) {
        let allFound = true;
        for (let d = 0; d < dayColumns.length; d++) {
          if (!visible[d] && dayColDataMasks[d].get(r))
            visible[d] = 1;
          if (!visible[d])
            allFound = false;
        }
        if (allFound) break;
      }
      for (let d = 0; d < dayColumns.length; d++) {
        const gridCol = grid.col(dayColumns[d].name);
        if (gridCol)
          gridCol.visible = visible[d] === 1;
      }
    });

    grid.onAfterDrawOverlay.subscribe(() => {
      if (!connectPointsInput.value)
        return;
      const selectedSubject = pivotedDf.getTag('selectedSubject');
      if (!selectedSubject)
        return;
      const enabledArms = getEnabledArms(pivotedDf);
      if (enabledArms.length === 0)
        return;

      const overlayCtx = grid.overlay.getContext('2d');
      if (!overlayCtx)
        return;

      const minCol = pivotedDf.col('~rowMin')!;
      const maxCol = pivotedDf.col('~rowMax')!;
      const visibleDayColNames = dayColNames.filter((n) => {
        const gc = grid.col(n);
        return gc && gc.visible;
      });
      if (visibleDayColNames.length < 2)
        return;

      const firstDayGridCol = grid.col(visibleDayColNames[0]);
      if (!firstDayGridCol)
        return;

      for (const visCell of firstDayGridCol.getVisibleCells()) {
        if (!visCell.isTableCell)
          continue;
        const gridRow = visCell.gridRow;
        const tableRowIdx = visCell.tableRowIndex;
        if (tableRowIdx == null)
          continue;

        const rowMin = minCol.get(tableRowIdx);
        const rowMax = maxCol.get(tableRowIdx);
        const range = rowMax - rowMin;
        const pts: {x: number, y: number}[] = [];

        for (const dayColName of visibleDayColNames) {
          const cell = grid.cell(dayColName, gridRow);
          if (!cell)
            continue;
          const value = cell.cell.value;
          if (!value)
            continue;

          let armData: any;
          try { armData = JSON.parse(value); }
          catch { continue; }

          const bounds = cell.bounds;
          const armWidth = (bounds.width - 2 * PADDING_X) / enabledArms.length;

          for (const arm of Object.keys(armData)) {
            const enabledIdx = enabledArms.indexOf(arm);
            if (enabledIdx < 0)
              continue;
            const entry = armData[arm];
            const subIdx = entry.subjects.indexOf(selectedSubject);
            if (subIdx < 0)
              continue;
            const px = bounds.x + PADDING_X + enabledIdx * armWidth + armWidth / 2;
            const v = entry.values[subIdx];
            const py = range === 0
              ? bounds.y + bounds.height / 2
              : bounds.y + bounds.height - PADDING_Y - ((v - rowMin) / range) * (bounds.height - 2 * PADDING_Y);
            pts.push({x: px, y: py});
            break;
          }
        }

        if (pts.length >= 2) {
          overlayCtx.strokeStyle = 'rgba(80,80,80,0.6)';
          overlayCtx.lineWidth = 1;
          overlayCtx.beginPath();
          overlayCtx.moveTo(pts[0].x, pts[0].y);
          for (let p = 1; p < pts.length; p++)
            overlayCtx.lineTo(pts[p].x, pts[p].y);
          overlayCtx.stroke();
        }
      }
    });

    restoreBrowsePanelOnRemoval();
  };

  return {df: pivotedDf, onTableViewAddedFunc: onTableViewAdded};
}