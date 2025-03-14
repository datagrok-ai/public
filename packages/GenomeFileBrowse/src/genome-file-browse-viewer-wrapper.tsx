import {
  createViewState,
  JBrowseLinearGenomeView,
} from "@jbrowse/react-linear-genome-view";

import React, { useEffect, useState } from "react";
import * as ReactDOM from "react-dom/client";

export function GenomeBrowseViewWrapper({ props, onUpdate }) {
  const [state, setState] = useState(createViewState(props));

  useEffect(() => {
    if (onUpdate) {
      onUpdate(setState);
    }
  }, [onUpdate]);

  return state?.session?.view ? (
    <JBrowseLinearGenomeView viewState={state} />
  ) : (
    <div>Loading...</div> // Prevents crash when session is undefined
  );
}

let component: { host: HTMLElement; loadData: (d: string) => void } | null =
  null;

// export function getGenomeFileBrowserComponent(host: HTMLElement, config: string) {
//   if (component)
//     return component;

//   // const data : any = JSON.parse(config);
//   // if(data.assemblies?.length > 0)
//   //   data.assembly = data.assemblies[0];

//   // const rElement = React.createElement(GenomeBrowseViewWrapper, data, null);
//   // const root = ReactDOM.createRoot(host);
//   // const loadDataObj = {loadData: (d: string) => {
//   //   const data : any = JSON.parse(config);
//   //   if(data.assemblies?.length > 0)
//   //     data.assembly = data.assemblies[0];
//   // }};
//   // root.render(rElement);
//   // component = {host: host, loadData: (d: string) => loadDataObj.loadData(d)};
//   // return component;

//   React.createElement(GenomeBrowseViewWrapper, {
//     props: { name: "GRCh38" },
//     onUpdate: (setState) => {
//       this.setGenomeState = setState; // Store setState function
//     },
//   })
// }

let setGenomeState = (p0: () => any) => {};
let genomeBiewComponent: any = undefined;

function initGenomeViewer(host: HTMLElement, config: string) {
  const root = ReactDOM.createRoot(host);

  const data: any = JSON.parse(config);
  if (data.assemblies?.length > 0) data.assembly = data.assemblies[0];

  root.render(<GenomeViewer data={data}/>);
  // root.render(
  //   React.createElement(GenomeBrowseViewWrapper, {
  //     props: data,
  //     onUpdate: (setState) => {
  //       setGenomeState = setState; // Store setState for external updates
  //     },
  //   })
  // );
}

function updateGenomeViewer(config: string) {
  if (typeof setGenomeState === "function") {
    const data: any = JSON.parse(config);
    if (data.assemblies?.length > 0) 
      data.assembly = data.assemblies[0];
    if (genomeBiewComponent) 
      genomeBiewComponent.setData(data);
  }
}

class GenomeViewer extends React.Component<{ data: any }, { data: any }> {
  constructor(props: { data: any }) {
    super(props); 
    this.state = {
      data: props.data,
    };
  }

  setData(data: any) {
    this.setState({ data });
  }

  render() {
    return <JBrowseLinearGenomeView viewState={this.state.data ?? {}} />;
  }
}

export function getGenomeViewer(host: HTMLElement, config: string) {
  if (component) return component;
  initGenomeViewer(host, config);
  component = { host: host, loadData: updateGenomeViewer };
  return component;
}
