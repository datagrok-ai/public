/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Editor } from "ketcher-react";
import { RemoteStructServiceProvider } from 'ketcher-core'
import React from 'react';
import ReactDOM from 'react-dom';

export let _package = new DG.Package();

const structServiceProvider = new RemoteStructServiceProvider(
  'foo_api_path'
  //process.env.REACT_APP_API_PATH!
)

//name: ketch
export function ketch() {


  let props = {
    //staticResourcesUrl: process.env.PUBLIC_URL,
    staticResourcesUrl: _package.webRoot,
    structServiceProvider: structServiceProvider
  };

  let host = ui.div([], { style: {height: '500px', width: '500px'}});
  // @ts-ignore
  let component = React.createElement(Editor, props, null);
  ReactDOM.render(component, host);

  ui.dialog()
    .add(host)
    .show();
}

// function createSketcher() {
//   return <Editor
//     staticResourcesUrl={_package.webRoot!}
//     structServiceProvider={structServiceProvider}/>;
// }
