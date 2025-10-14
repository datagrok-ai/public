import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getDataFrameDescription, getViewerDescriptionsString } from './utils';
import { askImpl } from './package';

type IViewerResponse = {
  viewerType: string;
  properties: {[index: string]: any};
}

/**
 * Examples of queries:
 * - plot width by height
 * - show me age distribution
 * - age distribution by disease
 * - plot age, height, and weight on the 3d scatter plot
 * - Plot width by height in a trellis plot split by race
 * - Plot weight by height for selected rows only
 * - form with all columns starting with "s"
 * */
export async function askAiTableView(view: DG.TableView, question: string) {
  const prompt = getPrompt(view, question);
  //console.log(prompt);
  const answer = await askImpl(prompt);
  //console.log(answer);
  //grok.shell.info(answer);

  const response: IViewerResponse = JSON.parse(answer.message!.content!);
  if (DG.Viewer.CORE_VIEWER_TYPES.includes(response.viewerType))
    view.addViewer(response.viewerType, response.properties);
}

export function askAiTableViewDialog(view: DG.TableView) {
  const input = ui.input.textArea('Question', {value: 'Plot width by height'});

  ui.dialog('Ask AI')
    .addInput('Question', input)
    .onOK(async () => {
      askAiTableView(view, input.value);
    })
    .show();
}

const getPrompt = (view: DG.TableView, question: string) => {
  return `Below is the structure of the dataframe a user is working with, and the interfaces for Datagrok visualizations. 
Your task is to identify which visualization to use, and which properties to set for what the user asks.
For column name properties, you can only use names from the dataframe.
Give the answer in the JSON format with the "viewerType" identifying viewer type, and
"properties" containing properties, like in this example: 
{"viewerType": "Histogram", "properties": {"valueColumnName": "Age", "markerSize": 15}}

User question: ${question}

DataFrame structure:
${JSON.stringify(getDataFrameDescription(view.dataFrame))}

Viewers:
${getViewerDescriptionsString()}`
}

