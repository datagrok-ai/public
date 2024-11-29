import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MmpInput} from './mmp-viewer';

export type MmpFilters = {
    activitySliderInputs: DG.InputBase[];
    activityValuesDivs: HTMLDivElement[];
    activityColorInputs: DG.InputBase[];
    activityActiveInputs: DG.InputBase[];
    filtersDiv: HTMLDivElement;
}

export function getMmpFilters(mmpInput: MmpInput, maxActs: number[]): MmpFilters {
  const activitySliderInputs = new Array<DG.InputBase>(maxActs.length);
  const activityValuesDivs = new Array<HTMLDivElement>(maxActs.length);
  const activityColorInputs = new Array<DG.InputBase>(maxActs.length);
  const activityActiveInputs = new Array<DG.InputBase>(maxActs.length);

  for (let i = 0; i < maxActs.length; i++) {
    const actName = mmpInput.activities.byIndex(i).name;
    const sliderInput = ui.input.slider(mmpInput.activities.byIndex(i).name,
      {value: maxActs[i] / 2, min: 0, max: maxActs[i]});
    const sliderInputValueDiv = ui.divText(sliderInput.stringValue, 'ui-input-description');
    sliderInput.addOptions(sliderInputValueDiv);
    sliderInput.root.classList.add('mmpa-slider-input');
    ui.tooltip.bind(sliderInput.captionLabel, `Select the cutoff by ${actName} difference`);
    ui.tooltip.bind(sliderInput.input, `${actName} value cutoff`);
    activitySliderInputs[i] = sliderInput;
    activityValuesDivs[i] = sliderInputValueDiv;
    const colorInput = ui.input.color('', {value: '#FF0000'});
    colorInput.root.classList.add('mmpa-color-input');
    activityColorInputs[i] = colorInput;
    const activeInput = ui.input.bool('', {value: true});
    activeInput.classList.add('mmpa-bool-input');
    activityActiveInputs[i] = activeInput;
  };


  const activitiesRoots: HTMLDivElement[] = new Array<HTMLDivElement>(activitySliderInputs.length);
  for (let i = 0; i < activitySliderInputs.length; i ++) {
    activitiesRoots[i] =
        ui.divH([activityActiveInputs[i].root, activityColorInputs[i].root, activitySliderInputs[i].root],
          {style: {paddingRight: '20px', height: '20px'}});
  }
  const activitiesDiv = ui.divV(activitiesRoots, 'mmpa-activity-filters');
  const filtersDiv = ui.divV([
    activitiesDiv,
  ], 'mmpa-slider-roots');

  return {activitySliderInputs, activityValuesDivs, activityColorInputs, activityActiveInputs,
    filtersDiv};
}
