import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  SeqAnnotation, SeqAnnotationHit, AnnotationCategory,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {getAnnotationColumnName, cacheAllRowAnnotations} from './annotations/annotation-manager';
import {_package} from '../package';

export function getRegionUI(col: DG.Column<string>): void {
  showGetRegionDialog(col);
}

/** Shows the Get Region dialog for the given column.
 *  When the user confirms, the region column is extracted, added to the dataframe, and
 *  {@link onRegionCreated} is called with it (if provided).
 *  Returns the dialog instance for further control. */
export function showGetRegionDialog(
  col: DG.Column<string>,
  onRegionCreated?: (regCol: DG.Column<string>) => void,
): DG.Dialog {
  const sh = _package.seqHelper.getSeqHandler(col);

  const nameInput = ui.input.string('Name', {value: ''});
  const startPositionInput = ui.input.choice('Start Position', {value: sh.posList[0], items: sh.posList,
    onValueChanged: () => updateNamePlaceholder()});
  const endPositionInput = ui.input.choice('End Position', {value: sh.posList[sh.posList.length - 1], items: sh.posList,
    onValueChanged: () => updateNamePlaceholder()});

  let selectedRegionName: string | null = null;

  const getDefaultName = (): string => {
    return selectedRegionName
      ? `${col.name}(${selectedRegionName}): ${startPositionInput.value}-${endPositionInput.value}`
      : `${col.name}:${startPositionInput.value}-${endPositionInput.value}`;
  };

  const updateNamePlaceholder = (): void => {
    if (!nameInput.value)
      nameInput.input.setAttribute('placeholder', getDefaultName());
  };
  updateNamePlaceholder();

  // Build region presets from annotations (new system) or .regions tag (legacy)
  const regionInput = _buildRegionPresetsInput(col, startPositionInput, endPositionInput, (regionName) => {
    selectedRegionName = regionName;
    updateNamePlaceholder();
  });

  const inputsList: DG.InputBase[] = [];
  if (regionInput) inputsList.push(regionInput);
  inputsList.push(nameInput, startPositionInput, endPositionInput);

  const dlg = ui.dialog({title: 'Get Region'}).add(ui.inputs(inputsList))
    .onOK(() => {
      const pi = DG.TaskBarProgressIndicator.create('Getting region...');
      try {
        const name: string = nameInput.value || getDefaultName();
        const regCol = getRegionDo(col, startPositionInput.value, endPositionInput.value, name);
        col.dataFrame.columns.add(regCol);
        regCol.setTag(DG.TAGS.CELL_RENDERER, 'sequence');
        onRegionCreated?.(regCol);
      } catch (err: any) {
        grok.shell.error(err.toString());
      } finally { pi.close(); }
    });
  dlg.show();
  return dlg;
}

/** Builds a Region preset dropdown from column annotations / legacy .regions tag.
 *  Returns null if the column has no annotated regions. */
function _buildRegionPresetsInput(
  col: DG.Column<string>,
  startInput: DG.InputBase<string | null>,
  endInput: DG.InputBase<string | null>,
  onRegionSelected?: (regionName: string | null) => void,
): DG.InputBase | null {
  type RegionPreset = { name: string, start: string, end: string };
  let regionList: RegionPreset[] | null = null;

  // New annotation system
  const annotationsTag: string | null = col.getTag(bioTAGS.annotations);
  if (annotationsTag) {
    try {
      const annotations = JSON.parse(annotationsTag);
      const structAnnots = annotations.filter(
        (a: any) => a.category === AnnotationCategory.Structure && a.start && a.end);
      if (structAnnots.length > 0) {
        regionList = structAnnots.map((a: any) => ({
          name: a.name, start: a.start, end: a.end,
        }));
      }
    } catch { /* ignore parse errors */ }
  }

  // Legacy .regions tag
  if (!regionList) {
    const regionsTagTxt: string | null = col.getTag(bioTAGS.regions);
    if (regionsTagTxt) {
      try { regionList = JSON.parse(regionsTagTxt); } catch { /* ignore */ }
    }
  }

  if (!regionList || regionList.length === 0) return null;

  const items = ['', ...regionList.map((r) => `${r.name}: ${r.start}-${r.end}`)];
  const regionInput = ui.input.choice('Region', {
    value: '', items: items,
    onValueChanged: (value: string) => {
      if (!value) {
        onRegionSelected?.(null);
        return;
      }
      const preset = regionList!.find((r) => `${r.name}: ${r.start}-${r.end}` === value);
      if (preset) {
        startInput.value = preset.start;
        endInput.value = preset.end;
        onRegionSelected?.(preset.name);
      }
    },
  });
  return regionInput;
}

/** {@link startPosName} and {@link endPosName} are according positionNames tag (or default ['1', '2',...]) */
export function getRegionDo(
  col: DG.Column<string>, startPosName: string | null, endPosName: string | null, name: string | null
): DG.Column<string> {
  const sh = _package.seqHelper.getSeqHandler(col);

  let startPosIdx: number | null = null;
  let endPosIdx: number | null = null;

  for (let posJ: number = 0; posJ < sh.posList.length; ++posJ) {
    if (sh.posList[posJ] == startPosName) startPosIdx = posJ;
    if (sh.posList[posJ] == endPosName) endPosIdx = posJ;
  }
  if (startPosIdx == null && startPosName !== null)
    throw new Error(`Start position ${startPosName} not found.`);
  if (endPosIdx == null && endPosName !== null)
    throw new Error(`End position ${endPosName} not found.`);

  if (sh.posList.length < endPosIdx!)
    throw new Error(`End position ${endPosIdx} exceeds positions length`);

  const regColName: string = !!name ? name : `${col.name}: (${startPosName ?? ''}-${endPosName ?? ''})`;

  // Try per-row extraction for unaligned data: find a matching structure annotation with per-row spans
  const perRowCol = _tryPerRowExtraction(col, sh, startPosName, endPosName, regColName);
  if (perRowCol) return perRowCol;

  // Fall back to column-level extraction (aligned/MSA data)
  const regCol = sh.getRegion(startPosIdx, endPosIdx, regColName);
  return regCol;
}

/** Attempts per-row region extraction using the companion annotation column.
 *  Returns null if no per-row data is available for the given start/end. */
function _tryPerRowExtraction(
  col: DG.Column<string>,
  sh: ReturnType<typeof _package.seqHelper.getSeqHandler>,
  startPosName: string | null,
  endPosName: string | null,
  regColName: string,
): DG.Column<string> | null {
  if (!startPosName || !endPosName || !col.dataFrame) return null;

  // Find matching structure annotation
  const annotTag = col.getTag(bioTAGS.annotations);
  if (!annotTag) return null;

  let annotations: SeqAnnotation[];
  try { annotations = JSON.parse(annotTag); } catch { return null; }

  const annot = annotations.find((a) =>
    a.category === AnnotationCategory.Structure && a.start === startPosName && a.end === endPosName);
  if (!annot) return null;

  // Check for companion annotation column with per-row region spans
  const annotColName = getAnnotationColumnName(col.name);
  let annotCol: DG.Column<string> | null = null;
  try { annotCol = col.dataFrame.columns.byName(annotColName) as DG.Column<string>; } catch { return null; }
  if (!annotCol) return null;

  const allRowData = cacheAllRowAnnotations(annotCol);
  const hasPerRowRegions = allRowData.some((rd) =>
    rd?.some((h: SeqAnnotationHit) => h.annotationId === annot.id && h.endPositionIndex != null));
  if (!hasPerRowRegions) return null;

  // Extract per-row using row-specific character indices
  const df = col.dataFrame;
  const regCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, regColName, df.rowCount);
  for (let i = 0; i < df.rowCount; i++) {
    const rowHits = allRowData[i];
    const regionHit = rowHits?.find((h: SeqAnnotationHit) =>
      h.annotationId === annot.id && h.endPositionIndex != null);
    if (regionHit) {
      const splitted = sh.getSplitted(i);
      const parts: string[] = [];
      for (let p = regionHit.positionIndex; p <= regionHit.endPositionIndex!; p++) {
        if (p < splitted.length)
          parts.push(splitted.getOriginal(p));
      }
      regCol.set(i, parts.join(sh.separator || ''));
    } else
      regCol.set(i, '');
  }
  return regCol;
}
