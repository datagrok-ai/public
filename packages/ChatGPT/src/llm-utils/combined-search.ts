/* eslint-disable max-len */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ModelType, OpenAIClient} from './openAI-client';
/**
 * The idea is to provide single AI search provider which will decide which LLM to use based on user prompt.
 * This class will not care about individual LLM implementations, nor the presence of api keys and other needed things.
 */
export class CombinedAISearchAssistant {
  private collectFunctions() {
    // these function should have descriptions, one input of a string prompt and an output of a widget. meta.useWhen must be defined as well
    const functions = DG.Func.find({meta: {role: 'aiSearchProvider'}}).filter((f) => f.inputs.length === 1 && f.inputs[0].type === 'string' &&
        f.description && f.outputs.length === 1 && (f.outputs[0].type as string) === 'widget' && f.options.useWhen);

    // create a map from names to functions
    const functionMap: {[name: string]: {func: DG.Func, friendlyName: string, inputName: string, order: number}} = {};
    for (const f of functions)
      functionMap[f.name] = {func: f, friendlyName: f.friendlyName ?? f.name, inputName: f.inputs[0].name, order: Number.parseInt(f.options.order ?? '0') ?? 0};
    return functionMap;
  }

  private functions: ReturnType<typeof CombinedAISearchAssistant.prototype.collectFunctions>;
  private constructor() {
    this.functions = this.collectFunctions();
  }
  private static _instance: CombinedAISearchAssistant | null = null;
  public static get instance(): CombinedAISearchAssistant {
    if (this._instance === null)
      this._instance = new CombinedAISearchAssistant();
    return this._instance;
  }

  public async getResultWidget(prompt: string): Promise<DG.Widget | null> {
    const functionOrder = await this.getRelevantFunctionOrder(prompt);
    console.log(`Function order for execution: ${functionOrder.join(', ')}`);
    const remainingFunctions = Object.keys(this.functions).filter((fn) => !functionOrder.includes(fn));
    let moreBtn: HTMLButtonElement | null = null;
    // create a lazy tabcontrol with all functions, but the chosen one first
    const tabcontrol = ui.tabControl();

    const addFunctionTab = (functionName: string, setAsCurrentPane: boolean = false) => {
      const funcInfo = this.functions[functionName];
      const pane = tabcontrol.addPane(funcInfo.friendlyName, () => {
        const wait = ui.wait(async () => {
          const inputParam: {[key: string]: any} = {};
          inputParam[funcInfo.inputName] = prompt;
          const result: DG.Widget | null = await funcInfo.func.apply(inputParam);
          if (result == null)
            return ui.divText('No result generated.');
          result.root.style.width = '100%';
          return result.root;
        });
        return wait;
      });
      if (setAsCurrentPane) tabcontrol.currentPane = pane;
      ui.tooltip.bind(pane.header, funcInfo.func.description);
    };

    const handleBtnClick = () => {
      const menu = ui.popupMenu();
      for (const functionName of remainingFunctions) {
        const funcFn = this.functions[functionName].friendlyName;
        menu.item(funcFn, () => {
          // Remove selected function from remaining list
          const index = remainingFunctions.indexOf(functionName);
          if (index !== -1) remainingFunctions.splice(index, 1);

          if (moreBtn && tabcontrol.header.contains(moreBtn))
            tabcontrol.header.removeChild(moreBtn);

          addFunctionTab(functionName, true);

          if (moreBtn && remainingFunctions.length > 0)
            tabcontrol.header.appendChild(moreBtn);
        });
      }
      menu.show();
    };

    for (const functionName of functionOrder)
      addFunctionTab(functionName);

    if (remainingFunctions.length) {
      moreBtn = ui.button('...', () => handleBtnClick(), 'Click to see other options');
      tabcontrol.header.appendChild(moreBtn);
    }

    tabcontrol.root.style.width = '100%';
    tabcontrol.root.style.height = 'unset';
    tabcontrol.root.style.minHeight = '300px';
    return DG.Widget.fromRoot(tabcontrol.root);
  }

  private lastUiprompt: string = '';
  public resetSearchUI() {
    this.lastUiprompt = '';
  }

  public async searchUI(prompt: string) {
    if (this.lastUiprompt === prompt)
      return;
    this.lastUiprompt = prompt;
    // add components to right places in ui
    const searchResultHost = document.getElementsByClassName('power-pack-search-host')[0];
    if (!searchResultHost) {
      console.warn('No search result host found for combined AI search.');
      return;
    }
    const spinner = ui.icons.spinner();
    spinner.style.color = 'var(--blue-1)';
    const searchingText = ui.divText('Grokking your question...', {style: {marginLeft: '8px'}});
    const loader = ui.divH([spinner, searchingText], {style: {alignItems: 'center', marginTop: '8px'}});
    const widgetDiv = ui.divV([loader], {style: {marginTop: '8px', marginBottom: '8px', alignItems: 'center', justifyContent: 'center', height: '100px'}});
    searchResultHost.prepend(widgetDiv);
    try {
      const widgetResult = await this.getResultWidget(prompt);
      if (!widgetResult)
        throw new Error('No result widget generated.');
      loader.remove();
      widgetDiv.style.removeProperty('align-items');
      widgetDiv.style.removeProperty('justify-content');
      widgetDiv.style.removeProperty('height');
      widgetDiv.appendChild(widgetResult.root);
    } catch (error: any) {
      console.error('Error during Combined AI search:', error);
      widgetDiv.appendChild(ui.divText(`Error during Combined AI search: ${error.message}`));
    } finally {
      loader.remove();
    }
  }

  private async getRelevantFunctionOrder(prompt: string) {
    const fullPrompt = `
        You are an intelligent assistant that helps to choose the best AI tool for a user's search query.
        Given the following available functions and their descriptions along with descriptions of when to use them, reply with JSON array of only relevant functions in correct order and nothing else.
        Make sure to reply WITH ONLY THE JSON PARSABLE ARRAY OF FUNCTION NAMES IN CORRECT ORDER (First being the best choice) and nothing else. no markdown, no backticks, no explanations nothing else.
        Leave out the functions only if they are super irrelevant, otherwise include them but make sure the best function is first in the list.
        Available functions:
        ${Object.values(this.functions).map((f) => `Function name: ${f.func.name}
        Description: ${f.func.description}
        Use when: ${f.func.options.useWhen}`).join('\n\n')}

        User query: ${prompt}

        `;
    const res = await OpenAIClient.getInstance().generalPromptCached(
      ModelType.Fast, 'You are a helpful assistant that selects the best function for a user query. reply with JSON array only and nothing else', fullPrompt);

    let order = Object.keys(this.functions);
    try {
      const parsed = JSON.parse(res) as string[];
      if (Array.isArray(parsed) && parsed.every((s) => typeof s === 'string' && this.functions[s]) && parsed.length > 0)
        order = parsed;
      else
        console.error('Invalid function order response:', res);
    } catch (error) {
      console.error('Error parsing function order response:', error, res);
    }
    return order;
  }
}
