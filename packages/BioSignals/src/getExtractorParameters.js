import * as ui from "datagrok-api/ui";

export function getExtractorParameters(extractorType) {
    switch (extractorType) {
        case 'Local energy':
            let winLen = ui.floatInput('Window Length', 2);
            winLen.setTooltip('Length of the window in seconds');
            let winStep = ui.floatInput('Window Step', 2);
            winStep.setTooltip('Shift of the window to start the next window');
            return {'win_len': winLen, 'win_step': winStep};
        case 'Beat from ECG':
            let bpm_max = ui.intInput('BPM max', 120);
            bpm_max.setTooltip('Maximal expected heart rate (in beats per minute), should be in range (1, 400]');
            let delta = ui.floatInput('Delta', 0);
            delta.setTooltip('Threshold for the peak detection (>=0). By default it is computed from the signal (adaptive thresholding)');
            let k = ui.floatInput('k', 0.7);
            k.setTooltip('Ratio (0,1) at which the signal range is multiplied (when delta = 0)');
            return {'bpm_max': bpm_max, 'delta': delta, 'k': k};
        case 'Phasic estimation':
            let delta1 = ui.floatInput('delta', 0.1);
            delta1.setTooltip('Minimum amplitude of the peaks in the driver');
            let t1 = ui.floatInput('t1', 0.75);
            t1.setTooltip('Value of the T1 parameter of the bateman function (>0)');
            let t2 = ui.floatInput('t2', 2);
            t2.setTooltip('Value of the T2 parameter of the bateman function (>0)');
            let grid_size = ui.floatInput('Grid size', 1);
            grid_size.setTooltip('Sampling size of the interpolation grid (>0)');
            let pre_max = ui.floatInput('Pre max', 2);
            pre_max.setTooltip('Duration (in seconds) of interval before the peak where to search the start of the peak (>0)');
            let post_max = ui.floatInput('Post max', 2);
            post_max.setTooltip('Duration (in seconds) of interval after the peak where to search the end of the peak (>0)');
            return {'delta': delta1, 't1': t1, 't2': t2, 'grid_size': grid_size, 'pre_max': pre_max, 'post_max': post_max};
        case 'BeatFromBP':
            let bpm_max1 = ui.intInput('BPM max', 120);
            bpm_max1.setTooltip('Maximal expected heart rate (in beats per minute), should be in range (1, 400]');
            let win_pre = ui.floatInput('win_pre', 0.25);
            win_pre.setTooltip('Portion (in seconds) to consider before the candidate beat position where to look for the beat, should be in range (0, 1]')
            let win_post = ui.floatInput('win_post', 0.05);
            win_post.setTooltip('Portion (in seconds) to consider after the candidate beat position where to look for the beat, should be in range (0, 1]');
            return {'bpm_max': bpm_max1, 'win_pre': win_pre, 'win_post': win_post};
    }
}