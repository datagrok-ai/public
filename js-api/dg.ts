import $ from 'cash-dom';
import * as dayjs from 'dayjs';
import * as wu from 'wu';

export * from './src/interfaces/d4';
export * from './src/const';
export * from './src/events';
export * from './src/dapi';
export * from './src/dataframe';
export * from './src/entities';
export * from './src/shell';
export * from './src/functions';
export * from './src/grid';
export * from './src/widgets';
export * from './src/views/view';
export * from './src/views/card_view';
export * from './src/views/files_view';
export * from './src/views/multi_view';
export * from './src/viewer';
export * from './src/docking';
export * from './src/wrappers_impl';
export * from './src/ui/wizard';
export {time, timeAsync, Utils, HtmlUtils, LruCache, Paint} from './src/utils';
export {ObjectHandler, EntityMetaDartProxy} from './ui';
export * from './src/data';
export * from './src/helpers';
export * from './src/logger';

export * from './src/chem';
export * from './src/ml';
export * from './src/utils';
export * from './src/ui/tree-view';

$(function () {
  (<any>window).$ = $;
  (<any>window).dayjs = dayjs;
  (<any>window).wu = wu;

  window.addEventListener("error", function (e) {
    e.preventDefault();
    e.stopPropagation();
    if (e.error.message == '[object ProgressEvent]')
      return;
    (<any>window).grok_Unhandled_Error(e.error.message ?? e.error, e.error.stack);

  });

});
