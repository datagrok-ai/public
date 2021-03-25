import * as ui from "datagrok-api/ui";

export function getIndicatorParameters(indicatorType) {
  switch (indicatorType) {
    case 'HRV time domain':
      let step = ui.floatInput('Step', 5);
      step.setTooltip('Time distance between subsequent segments');
      let width = ui.floatInput('Window width', 10);
      width.setTooltip('Time window width');
      return {'step': step, 'width': width};
    case 'HRV frequency domain':
      let resampling_frequency = ui.floatInput('Resample to', 4);
      resampling_frequency.setTooltip('Resample RR intervals to be able to apply FFT (it requires constant time differences between samples)');
      let step2 = ui.floatInput('Step', 5);
      step2.setTooltip('Time distance between subsequent segments');
      let width2 = ui.floatInput('Window width', 10);
      width2.setTooltip('Time window width');
      return {'resampling_frequency': resampling_frequency, 'step': step2, 'width': width2};
    case 'HRV nonlinear domain':
      let step1 = ui.floatInput('Step', 5);
      step1.setTooltip('Time distance between subsequent segments');
      let width1 = ui.floatInput('Window width', 10);
      width1.setTooltip('Time window width');
      return {'step': step1, 'width': width1};
  }
}