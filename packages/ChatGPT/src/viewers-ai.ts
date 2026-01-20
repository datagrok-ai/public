import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getCurrentViewersString, getDataFrameDescription, getViewerDescriptionsString} from './utils';
import {ModelType, OpenAIClient} from './llm-utils/openAI-client';

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
  // TODO: substitute with the OpenAI client call
  const res = OpenAIClient.getInstance();
  const answer = await res.generalPromptCached(ModelType.Fast, '', prompt);
  //console.log(answer);
  //grok.shell.info(answer);

  const response: IViewerResponse = JSON.parse(answer);
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
  return `
  Below is the structure of the dataframe a user is working with, and the interfaces for Datagrok visualizations. 
Your task is to identify which visualization to use, and which properties to set for what the user asks.
For column name properties, you can only use names from the dataframe.
Give the answer in the JSON format that can be directly parsed (no annotations or anything like that)
 with the "viewerType" identifying viewer type, and
"properties" containing properties, like in this example: 
{"viewerType": "Histogram", "properties": {"valueColumnName": "Age", "markerSize": 15}}

${getCurrentViewersString(view)}

User question: ${question}

DataFrame structure:
${JSON.stringify(getDataFrameDescription(view.dataFrame))}

Viewers:
${getViewerDescriptionsString()}`;
};

