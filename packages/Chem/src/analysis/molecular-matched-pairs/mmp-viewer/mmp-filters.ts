import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type MmpFilters = {
    activitySliderInputs: DG.InputBase[];
    activityValuesDivs: HTMLDivElement[];
    activityColorInputs: DG.InputBase[];
    activityActiveInputs: DG.InputBase[];
    filtersDiv: HTMLDivElement;
}

export function getMmpFilters(activities: string[], maxActs: number[]): MmpFilters {
  const activitySliderInputs = new Array<DG.InputBase>(maxActs.length);
  const activityValuesDivs = new Array<HTMLDivElement>(maxActs.length);
  const activityColorInputs = new Array<DG.InputBase>(maxActs.length);
  const activityActiveInputs = new Array<DG.InputBase>(maxActs.length);

  for (let i = 0; i < maxActs.length; i++) {
    const actName = activities[i];
    //set activities filters to values a bit higher than 0.5 from max value: maxValue/1.7
    const sliderInput = ui.input.slider(activities[i],
      {value: maxActs[i] / 1.7, min: 0, max: maxActs[i]});
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
