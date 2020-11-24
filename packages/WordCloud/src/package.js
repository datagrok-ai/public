/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import { WordCloudViewer } from './word-cloud-viewer.js';
import '../css/word-cloud.css';

export const _package = new DG.Package();

//name: Word Cloud
//description: Creates a word cloud
//tags: viewer
//output: viewer result
export function wordcloud() {
    return new WordCloudViewer();
}
