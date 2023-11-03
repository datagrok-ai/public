/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ExternalPluginUI} from '../view/app-ui';
import {ColoredTextInput} from '../view/utils/colored-input/colored-text-input';
import {highlightInvalidSubsequence} from '../view/utils/colored-input/input-painters';
import {codesToSymbolsDictionary} from '../model/data-loading-utils/json-loader';
import {MERMADE} from './const';

export async function getExternalAppViewFactories(): Promise<{[name: string]: () => DG.View} | undefined> {
  const externalPluginData = {
  [MERMADE.FUNCTION_NAME]: {
      tabName: MERMADE.TAB_NAME,
      parameters: getMerMadeParameters()
    },
  }

  const result: {[tabName: string]: () => DG.View} = {};

  for (const [pluginName, data] of Object.entries(externalPluginData)) {
    let div: HTMLDivElement;
    try {
      div = await grok.functions.call(pluginName, data.parameters);
    } catch (err) {
      console.error(`Plugin ${pluginName} not loaded, error:`, err)
      div = ui.divText('error loading');
    }
    const pluginUI = new ExternalPluginUI(data.tabName, div);

    // intentonally don't await for the promise
    pluginUI.initView();

    result[data.tabName] = () => pluginUI.getView();
  }
  return result;
}

function getMerMadeParameters(): {[name: string]: any} {
  const base = ui.textInput('', '');
  const input = new ColoredTextInput(base, highlightInvalidSubsequence);

  return {
    coloredInput: input,
    codes: codesToSymbolsDictionary
  }
}
