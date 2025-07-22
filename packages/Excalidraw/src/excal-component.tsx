import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Excalidraw} from '@excalidraw/excalidraw';
import React, { useState } from "react";
// import ReactDOM from "";
import * as ReactDOM from 'react-dom/client';

function App(props: {obj: any}) {
    const [excalidrawAPI, setExcalidrawAPI] = useState<any>(null);
    props.obj ??= {};
    props.obj.get = () => excalidrawAPI; 
    return (
      <>
        <div style={{ height: "100%" }}>
          <Excalidraw excalidrawAPI={(api) => setExcalidrawAPI(api)}/>
        </div>
      </>
    );
  }

export function renderExcalidraw() {
    const actualRoot = ui.div([], 'd4-excalidraw-root');
    const ketcherHost = ui.div([], 'd4-excalidraw-host');
    actualRoot.style.height = '100%';
    ketcherHost.style.height = 'calc(100% - 20px)';
    const obj: any = {get: () => null};
    const component = React.createElement(App, {obj}, null);
    const root = ReactDOM.createRoot(ketcherHost);
    root.render(component);
    actualRoot.appendChild(ketcherHost);
    return {actualRoot, getApi: () => obj.get()};
}