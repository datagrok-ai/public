/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Editor } from "ketcher-react";
import { RemoteStructServiceProvider } from 'ketcher-core'
import React from 'react';
import ReactDOM from 'react-dom';

export let _package = new DG.Package();

function react(reactComponent: React.DOMElement<any, any> | Array<React.DOMElement<any, any>> | React.CElement<any, any> | Array<React.CElement<any, any>> | React.ReactElement | Array<React.ReactElement>): DG.Widget {
  let widget = DG.Widget.fromRoot(ui.div());
  // @ts-ignore
  ReactDOM.render(reactComponent, widget.root);
  return widget;
}

//name: ketch
export function ketch() {

  const structServiceProvider = new RemoteStructServiceProvider(
    'foo_api_path'
    //process.env.REACT_APP_API_PATH!
  )

  let props = {
    //staticResourcesUrl: process.env.PUBLIC_URL,
    staticResourcesUrl: _package.webRoot,
    structServiceProvider: structServiceProvider
  };

  let sketcherWidget = react(
    // @ts-ignore
    React.createElement(Editor, props, null)
  );

  ui.dialog()
    .add(sketcherWidget)
    .show();
}
