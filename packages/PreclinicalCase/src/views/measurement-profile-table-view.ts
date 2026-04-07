import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {createAllMeasurementsDf, createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SEX, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';
import {createSubjectProfileView} from './subject-profile-view';
import {StudyTableViewParams} from '../utils/views-creation-utils';
import {restoreBrowsePanelOnRemoval} from '../utils/utils';

function parseIso8601DurationToDays(duration: string): number | null {
  const m = duration.match(/^P(?:(\d+)W)?(?:(\d+)D)?$/);
  if (!m)
    return null;
  return (m[1] ? parseInt(m[1]) * 7 : 0) + (m[2] ? parseInt(m[2]) : 0);
}

function getRecoveryStartDay(studyId: string): number | null {
  const study = studies[studyId];
  const dosingStr = study.config.other?.['Dosing Duration'];
  if (dosingStr) {
    const dosingDays = parseIso8601DurationToDays(dosingStr);
    if (dosingDays !== null)
      return dosingDays + 1;
  }
  const te = study.domains.te;
  if (te) {
    const elementCol = te.col('ELEMENT');
    const testrlCol = te.col('TESTRL');
    if (elementCol && testrlCol) {
      for (let i = 0; i < te.rowCount; i++) {
        const element = elementCol.get(i);
        if (element && element.toLowerCase().includes('recovery')) {
          const testrl = testrlCol.get(i);
          if (testrl) {
            const days = parseIso8601DurationToDays(testrl);
            if (days !== null)
              return days;
          }
        }
      }
    }
  }
  return null;
}

export function createMeasurementProfileTableView(studyId: string): StudyTableViewParams {
  let resDf = createAllMeasurementsDf(studyId);

  if (!resDf)
    return {df: DG.DataFrame.create()};

  const columnsFromDm = studies[studyId].domains.dm!.columns.names().filter((it) => it !== SUBJECT_ID);
  grok.data.joinTables(resDf, studies[studyId].domains.dm!, [SUBJECT_ID], [SUBJECT_ID], null,
    columnsFromDm, DG.JOIN_TYPE.LEFT, true);

  const intVisitDay = (resDf as DG.DataFrame).col('visit_day')!.convertTo( DG.TYPE.INT);
  (resDf as DG.DataFrame).columns.remove('visit_day');
  (resDf as DG.DataFrame).columns.add(intVisitDay);

  (resDf as DG.DataFrame).name = 'Measurements';
  (resDf as DG.DataFrame).rows.filter((row) => row.test == (resDf as DG.DataFrame).get('test', 0));

  const recoveryStartDay = getRecoveryStartDay(studyId);
  if (recoveryStartDay !== null) {
    resDf.meta.formulaLines.addLine({
      formula: `\${visit_day} = ${recoveryStartDay}`,
      color: '#C83C3C',
      width: 1.5,
      style: 'dashed',
      title: 'Recovery start',
    });
  }

  const armCol = resDf.col(PLANNED_TRT_ARM);
  const subjectCol = resDf.col(SUBJECT_ID);
  const armColors: {[arm: string]: string} = {};
  if (armCol) {
    for (const cat of armCol.categories)
      armColors[cat] = DG.Color.toHtml(DG.Color.getCategoryColor(armCol, cat));
  }
  const subjectArmColors: {[subject: string]: string} = {};
  if (armCol && subjectCol) {
    for (let i = 0; i < resDf.rowCount; i++) {
      const subject = subjectCol.get(i);
      const arm = armCol.get(i);
      if (subject && arm && !subjectArmColors[subject])
        subjectArmColors[subject] = armColors[arm] ?? '#888888';
    }
  }

  const armCodeToName = studies[studyId].armCodeToName;
  const getArmCode = (arm: string) => Object.keys(armCodeToName).find((c) => armCodeToName[c] === arm) ?? arm;

  const armBitSets: {[arm: string]: DG.BitSet} = {};
  if (armCol) {
    for (const cat of armCol.categories)
      armBitSets[cat] = DG.BitSet.create(resDf.rowCount, (i) => armCol.get(i) === cat);
  }
  let baseFilter = resDf.filter.clone();
  const disabledArms = new Set<string>();
  let isArmFiltering = false;

  let trellisPlot: DG.Viewer | null = null;
  const onTableViewAdded = async (tableView: DG.TableView) => {
    const colorBy = ui.input.choice('Color by', {
      items: ['Subject', 'Treatment arm'],
      value: 'Subject',
      onValueChanged: () => {
        const subjCol = tableView.dataFrame.col(SUBJECT_ID);
        if (!subjCol) return;
        if (colorBy.value === 'Treatment arm') {
          subjCol.meta.colors.setCategorical(subjectArmColors);
          legendDiv.style.display = '';
        } else {
          subjCol.meta.colors.setCategorical({});
          legendDiv.style.display = 'none';
          if (disabledArms.size > 0)
            resetArmFilter();
        }
        tableView.dataFrame.fireValuesChanged();
      },
    });
    colorBy.root.style.display = 'none';

    const splitBy = ui.input.choice('Split by', {
      items: ['Treatment arm', 'Subject'],
      value: 'Treatment arm',
      onValueChanged: () => {
        const isSplitBySubject = splitBy.value === 'Subject';
        trellisPlot?.setOptions({
          innerViewerLook: {
            splitColumnNames: isSplitBySubject ? [SUBJECT_ID] : [PLANNED_TRT_ARM],
          },
        });
        colorBy.root.style.display = isSplitBySubject ? '' : 'none';
        if (!isSplitBySubject) {
          const subjCol = tableView.dataFrame.col(SUBJECT_ID);
          if (subjCol) {
            subjCol.setTag(DG.TAGS.COLOR_CODING_TYPE, '');
            subjCol.setTag(DG.TAGS.COLOR_CODING_CATEGORICAL, '');
          }
          legendDiv.style.display = 'none';
          if (disabledArms.size > 0)
            resetArmFilter();
          colorBy.value = 'Subject';
        }
      },
    });

    const armLabels: {arm: string, el: HTMLElement, color: string}[] = [];
    const legendDiv = ui.div([], {classes: 'preclinical-case-legend'});
    legendDiv.style.display = 'none';

    function updateArmFilter() {
      isArmFiltering = true;
      const df = resDf as DG.DataFrame;
      df.filter.copyFrom(baseFilter, false);
      if (disabledArms.size > 0 && armCol) {
        const enabledMask = DG.BitSet.create(df.rowCount);
        for (const cat of armCol.categories) {
          if (!disabledArms.has(cat))
            enabledMask.or(armBitSets[cat], false);
        }
        df.filter.and(enabledMask);
      } else {
        df.filter.fireChanged();
      }
      isArmFiltering = false;
    }

    function resetArmFilter() {
      disabledArms.clear();
      for (const entry of armLabels) {
        entry.el.style.color = entry.color;
        entry.el.style.opacity = '';
      }
      updateArmFilter();
    }

    if (armCol) {
      for (const arm of armCol.categories) {
        const color = armColors[arm] ?? '#888888';
        const armCode = getArmCode(arm);
        const label = ui.div([armCode], {classes: 'preclinical-case-legend-label', style: {color}});
        label.style.cursor = 'pointer';
        ui.tooltip.bind(label, () => {
          const tip = ui.divText(arm);
          tip.style.color = disabledArms.has(arm) ? '#aaa' : color;
          return tip;
        });
        label.addEventListener('click', () => {
          if (disabledArms.has(arm)) {
            disabledArms.delete(arm);
            label.style.color = color;
            label.style.opacity = '';
          } else {
            disabledArms.add(arm);
            label.style.color = '#ccc';
            label.style.opacity = '0.5';
          }
          updateArmFilter();
        });
        legendDiv.appendChild(label);
        armLabels.push({arm, el: label, color});
      }
    }

    resDf.onFilterChanged.subscribe(() => {
      if (isArmFiltering) return;
      baseFilter = resDf.filter.clone();
      if (disabledArms.size > 0)
        updateArmFilter();
    });

    const hasChange = tableView.dataFrame.col('change') !== null;
    const hasPctChange = tableView.dataFrame.col('pct_change') !== null;
    const valueItems = ['Absolute'];
    if (hasChange)
      valueItems.push('Change from baseline');
    if (hasPctChange)
      valueItems.push('% change from baseline');
    const valueType = ui.input.choice('Value', {
      items: valueItems,
      value: 'Absolute',
      onValueChanged: () => {
        const yCol = valueType.value === '% change from baseline' ? 'pct_change' :
          valueType.value === 'Change from baseline' ? 'change' : 'result';
        trellisPlot?.setOptions({
          innerViewerLook: {yColumnNames: [yCol]},
        });
      },
    });

    const showRecoveryLineInput = recoveryStartDay !== null
      ? ui.input.bool('Recovery line', {value: true, onValueChanged: (value) => {
          const lines = resDf.meta.formulaLines.items;
          const idx = lines.findIndex((l) => l.title === 'Recovery start');
          if (idx >= 0)
            resDf.meta.formulaLines.updateAt(idx, {...lines[idx], visible: value});
        }})
      : null;

    const ribbonItems = [valueType.root, splitBy.root, colorBy.root, legendDiv];
    if (showRecoveryLineInput)
      ribbonItems.push(showRecoveryLineInput.root);
    tableView.setRibbonPanels([ribbonItems]);

    trellisPlot = await DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, tableView.dataFrame, {
      xColumnNames: [SEX],
      yColumnNames: ['test'],
      viewerType: DG.VIEWER.LINE_CHART,
      innerViewerLook: {
        showXAxis: true,
        aggrType: DG.STATS.AVG,
        splitColumnNames: [PLANNED_TRT_ARM],
        showDataframeFormulaLines: true,
      },
    });

    let subjectView: DG.DockNode | null = null;
    tableView.dataFrame.onCurrentRowChanged.subscribe((_) => {
      if (tableView.dataFrame.currentRowIdx === -1)
        return;
      const look = trellisPlot?.getOptions().look;
      if (look?.viewerType === DG.VIEWER.LINE_CHART && look.innerViewerLook.splitColumnNames.includes(SUBJECT_ID)) {
        const v = createSubjectProfileView(studyId, tableView.dataFrame,
          tableView.dataFrame.get(SUBJECT_ID, tableView.dataFrame.currentRowIdx));

        if (subjectView)
          (tableView as DG.TableView).dockManager.close(subjectView);
        subjectView = (tableView as DG.TableView).dockManager.dock(v, DG.DOCK_TYPE.DOWN);
      }
    });

    tableView.dockManager.dock(trellisPlot, DG.DOCK_TYPE.FILL);

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL,
      column: 'test',
      columnName: 'test',
    });
    restoreBrowsePanelOnRemoval();
  };

  return {df: resDf, onTableViewAddedFunc: onTableViewAdded};
}
