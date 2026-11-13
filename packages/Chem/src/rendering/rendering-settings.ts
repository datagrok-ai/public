/**
 * The "Rendering" package settings (package.json `properties`): drawing defaults of every molecule Chem draws, in both
 * engines, which a molecule column overrides for its cells by a tag (Chemistry | Rendering).
 */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chemSettings} from '../engine/chem-engine';
import {hexToPercentRgb} from '../utils/chem-common';

export interface RenderingSetting {
  /** The package property, whose `friendlyName` in package.json is {@link caption}. */
  property: string;
  /** The molecule column's tag that overrides it. */
  tag: string;
  /** The value when the package has none (a package version without the property): package.json's default. */
  defaultValue: string;
  caption: string;
  description: string;
}

export const HIGHLIGHT_COLOR: RenderingSetting = {
  property: 'HighlightColor', tag: '.%chem-highlight-color', defaultValue: '#ffb3b3', caption: 'Highlight color',
  description: 'Color of highlighted substructures that have no color of their own, such as the matches of the ' +
    'substructure filter',
};
export const UNSPECIFIED_STEREO_IS_UNKNOWN: RenderingSetting = {
  property: 'UnspecifiedStereoIsUnknown', tag: '.%chem-unspecified-stereo-is-unknown', defaultValue: 'true',
  caption: 'Mark unspecified stereo',
  description: 'Draw stereocenters with unspecified configuration with a wavy bond, and double bonds with ' +
    'unspecified E/Z configuration as crossed bonds',
};
export const SIMPLIFIED_STEREO_GROUP_LABEL: RenderingSetting = {
  property: 'SimplifiedStereoGroupLabel', tag: '.%chem-simplified-stereo-group-label', defaultValue: 'true',
  caption: 'Simplified stereo label',
  description: 'When all stereocenters of a molecule are in one AND or OR enhanced stereo group, label the ' +
    'molecule once ("AND enantiomer", "OR enantiomer") instead of each stereocenter ("and1", "or1")',
};

/** The package's value of `setting`, as a tag holds it. */
export function renderingDefault(setting: RenderingSetting): string {
  return String(chemSettings()[setting.property] ?? setting.defaultValue);
}

/** `col`'s value of `setting`: its tag, else the package's. */
export function renderingValue(setting: RenderingSetting, col?: DG.Column | null): string {
  return col?.getTag(setting.tag) || renderingDefault(setting);
}

/** Sets `col`'s override of `setting`; the package's value removes it, so the column follows the package again. */
export function setRenderingValue(setting: RenderingSetting, col: DG.Column, value: string): void {
  if (value === renderingDefault(setting))
    col.tags.delete(setting.tag);
  else
    col.setTag(setting.tag, value);
  col.dataFrame?.fireValuesChanged();
}

/** RDKit's draw options of the Rendering settings for `col`'s cells, or the package's without a column. */
export function renderingOptions(col?: DG.Column | null): {[key: string]: any} {
  return {
    highlightColour: hexToPercentRgb(renderingValue(HIGHLIGHT_COLOR, col), 1) ??
      hexToPercentRgb(HIGHLIGHT_COLOR.defaultValue, 1),
    unspecifiedStereoIsUnknown: renderingValue(UNSPECIFIED_STEREO_IS_UNKNOWN, col) === 'true',
    simplifiedStereoGroupLabel: renderingValue(SIMPLIFIED_STEREO_GROUP_LABEL, col) === 'true',
  };
}

/** The inputs of `col`'s overrides (Chemistry | Rendering). */
export function renderingInputs(col: DG.Column): DG.InputBase[] {
  const color = ui.input.color(HIGHLIGHT_COLOR.caption, {value: renderingValue(HIGHLIGHT_COLOR, col),
    tooltipText: HIGHLIGHT_COLOR.description,
    onValueChanged: (value) => setRenderingValue(HIGHLIGHT_COLOR, col, value)});
  const flag = (setting: RenderingSetting) => ui.input.bool(setting.caption, {
    value: renderingValue(setting, col) === 'true', tooltipText: setting.description,
    onValueChanged: (value) => setRenderingValue(setting, col, value.toString())});
  return [color, flag(UNSPECIFIED_STEREO_IS_UNKNOWN), flag(SIMPLIFIED_STEREO_GROUP_LABEL)];
}
