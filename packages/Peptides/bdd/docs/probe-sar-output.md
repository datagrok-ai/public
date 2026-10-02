# probe-sar.mjs output (2026-09-11, stand localhost:8888, peptides.csv, Launch SAR from the Peptides pane)

```
0.4s table open
{
  rows: 647,
  cols: [ 'ID', 'AlignedSequence', 'IC50' ],
  semType: 'Macromolecule',
  units: 'separator',
  renderer: 'sequence'
}
0.5s Peptides pane
panes: pane-Details | pane-Filter | pane-Actions | pane-Colors | pane-Style | pane-Settings | pane-Plots | pane-Advanced | pane-Peptides | pane-Bioinformatics
4.2s pane expanded
pane names: div[div-section--Peptides] rect[pan-handle] circle[min-handle] circle[max-handle] div[input-host-Activity] div[input-Activity] i[icon-plus] i[icon-edit] div[input-host-Scaling] select[input-Scaling] div[input-host-Clusters] div[input-Clusters] i[icon-plus] i[icon-edit] div[input-host-Generate-clusters] input[input-Generate-clusters] button[button-Launch-SAR] div[viewer-Histogram] canvas[canvas] div[div-column-combobox-value] i[icon-plus] i[icon-edit] div[div-column-combobox-split] i[icon-plus] i[icon-edit] input[input-Stack-split-categories] input[input-Filter-out-missing-values] div[input-host-Similarity-Threshold] input[input-Similarity-Threshold] div[input-host-Inflation-Factor] input[input-Inflation-Factor] div[input-host-Max-Iterations] input[input-Max-Iterations] div[input-host-Min-Cluster-Size] input[input-Min-Cluster-Size] div[input-host-Distance-Function] select[input-Distance-Function] div[input-host-Fingerprint-Type] select[input-Fingerprint-Type] div[input-host-Gap-Open-Penalty] input[input-Gap-Open-Penalty] div[input-host-Gap-Extend-Penalty] input[input-Gap-Extend-Penalty] div[input-host-Use-WebGPU] input[input-Use-WebGPU]
pane weblogo: {
  host: true,
  rootName: null,
  widget: 'L',
  hasStatus: true,
  canvases: 2
}
4.7s model
38.2s 4 viewers
{
 "view": "Table",
 "viewers": [
  {
   "type": "Grid",
   "ctor": "Grid",
   "rootName": "viewer-Grid",
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "Grid"
  },
  {
   "type": "Sequence Variability Map",
   "ctor": "I",
   "rootName": null,
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "I"
  },
  {
   "type": "Most Potent Residues",
   "ctor": "O",
   "rootName": null,
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "O"
  },
  {
   "type": "MCL",
   "ctor": "mn",
   "rootName": null,
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "mn"
  },
  {
   "type": "Logo Summary Table",
   "ctor": "T",
   "rootName": null,
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "T"
  },
  {
   "type": "Grid",
   "ctor": "Grid",
   "rootName": "viewer-Grid",
   "status": "function",
   "pending": "boolean",
   "onRendered": "undefined",
   "found": "Grid"
  }
 ],
 "roots": [
  "viewer-Grid 225x471 in=-",
  "viewer-MCL 225x447 in=-",
  "viewer-Scatter-plot 225x447 in=viewer-MCL",
  "viewer-Sequence-Variability-Map 315x448 in=-",
  "viewer-Grid 315x425 in=viewer-Sequence-Variability-Map",
  "viewer-Most-Potent-Residues 135x448 in=-",
  "viewer-Grid 135x448 in=viewer-Most-Potent-Residues",
  "viewer-Logo-Summary-Table 453x540 in=-",
  "viewer-Grid 453x540 in=viewer-Logo-Summary-Table",
  "viewer-Grid 453x355 in=-"
 ],
 "ribbon": [
  "div name=null aria=null cls=d4-ribbon-item",
  "div name=view selector aria=null cls=d4-combo-popup",
  "div name=Home-host aria=null cls=",
  "div name=Home aria=null cls=d4-icon-text-small d4-list-item",
  "i name=null aria=null cls=grok-icon",
  "i name=icon-times aria=Close view cls=grok-icon fal fa-times",
  "div name=Table-host aria=null cls=",
  "div name=Table aria=null cls=d4-icon-text-small d4-list-item",
  "i name=icon-view-layout aria=null cls=grok-icon svg-icon svg-view-layout",
  "i name=icon-times aria=Close view cls=grok-icon fal fa-times",
  "i name=null aria=null cls=grok-icon",
  "i name=icon-window-maximize aria=null cls=grok-icon fal fa-window-maximize",
  "div name=div-view-name aria=null cls=grok-ns-entity-name d4-ribbon-name",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-plus aria=Add table to workspace cls=grok-icon fal fa-plus",
  "div name=null aria=null cls=d4-ribbon-item no-hover",
  "button name=button-Save aria=null cls=ui-btn ui-btn-ok ui-btn-raised",
  "i name=icon-cloud-upload aria=null cls=grok-icon fal fa-cloud-upload",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-arrow-to-bottom aria=null cls=grok-icon fal fa-arrow-to-bottom",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=null aria=Add viewer cls=grok-icon svg-icon svg-add-viewer",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-filter aria=Toggle filters cls=grok-icon far fa-filter grok-icon-filter",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-select-all aria=All cls=grok-icon svg-icon svg-select-all",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-select-none aria=None cls=grok-icon svg-icon svg-select-none",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-invert-selection aria=Invert cls=grok-icon svg-icon svg-invert-selection",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-remove-selected-rows aria=null cls=grok-icon svg-icon svg-remove-selected-rows d4-disabled",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-remove-selected-columns aria=null cls=grok-icon svg-icon svg-remove-selected-columns d4-disabled",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=icon-add-new-column aria=Add New Column... cls=grok-icon svg-icon svg-add-new-column",
  "div name=null aria=null cls=d4-ribbon-item",
  "i name=null aria=Peptides analysis settings cls=grok-icon fal fa-wrench",
  "i name=icon-font-icon-menu aria=null cls=grok-icon grok-font-icon-menu",
  "i name=icon-times aria=null cls=grok-icon fal fa-times"
 ],
 "grids": [
  {
   "name": "viewer-Grid",
   "widget": "Grid",
   "type": "Grid",
   "areas": [
    "x scroll slider",
    "x scroll handle",
    "y scroll slider",
    "y scroll handle",
    "cell 1 of 10",
    "cell 2 of 10"
   ],
   "values": [
    "rows shown",
    "rows",
    "pinned rows",
    "sort column",
    "sort direction",
    "is heatmap",
    "heatmap colors",
    "global color scaling",
    "max heatmap columns",
    "row height",
    "col labels orientation",
    "effective col labels orientation"
   ]
  },
  {
   "name": "viewer-Grid",
   "inside": "viewer-Sequence-Variability-Map",
   "widget": "Grid",
   "type": "Grid",
   "areas": [
    "x scroll slider",
    "x scroll handle",
    "y scroll slider",
    "y scroll handle",
    "cell 2 of AAR",
    "cell 3 of AAR"
   ],
   "values": [
    "rows shown",
    "rows",
    "pinned rows",
    "sort column",
    "sort direction",
    "is heatmap",
    "heatmap colors",
    "global color scaling",
    "max heatmap columns",
    "row height",
    "col labels orientation",
    "effective col labels orientation"
   ]
  },
  {
   "name": "viewer-Grid",
   "inside": "viewer-Most-Potent-Residues",
   "widget": "Grid",
   "type": "Grid",
   "areas": [
    "x scroll slider",
    "x scroll handle",
    "cell 1 of Pos",
    "cell 2 of Pos",
    "cell 3 of Pos",
    "cell 4 of Pos"
   ],
   "values": [
    "rows shown",
    "rows",
    "pinned rows",
    "sort column",
    "sort direction",
    "is heatmap",
    "heatmap colors",
    "global color scaling",
    "max heatmap columns",
    "row height",
    "col labels orientation",
    "effective col labels orientation"
   ]
  },
  {
   "name": "viewer-Grid",
   "inside": "viewer-Logo-Summary-Table",
   "widget": "Grid",
   "type": "Grid",
   "areas": [
    "x scroll slider",
    "x scroll handle",
    "cell 2 of Cluster",
    "header Cluster",
    "column resizer Cluster",
    "cell 2 of Members"
   ],
   "values": [
    "rows shown",
    "rows",
    "pinned rows",
    "sort column",
    "sort direction",
    "is heatmap",
    "heatmap colors",
    "global color scaling",
    "max heatmap columns",
    "row height",
    "col labels orientation",
    "effective col labels orientation"
   ]
  },
  {
   "name": "viewer-Grid",
   "widget": "Grid",
   "type": "Grid",
   "areas": [
    "x scroll slider",
    "x scroll handle",
    "header Activity",
    "column resizer Activity",
    "header AlignedSequence",
    "column resizer AlignedSequence"
   ],
   "values": [
    "rows shown",
    "rows",
    "pinned rows",
    "sort column",
    "sort direction",
    "is heatmap",
    "heatmap colors",
    "global color scaling",
    "max heatmap columns",
    "row height",
    "col labels orientation",
    "effective col labels orientation"
   ]
  }
 ]
}
svm inputs: div[input-host-Mutation-Cliffs] input[input-Mutation-Cliffs] div[input-host-Invariant-Map] input[input-Invariant-Map] input[input-Search] div[viewer-Grid] svg[x-slider] rect[pan-handle] circle[min-handle] circle[max-handle] svg[y-slider] rect[pan-handle] circle[min-handle] circle[max-handle] canvas[canvas] canvas[overlay] i[icon-font-icon-settings] i[icon-plus] i[icon-minus] i[icon-font-icon-menu] div[viewer-Grid] svg[x-slider] rect[pan-handle] circle[min-handle] circle[max-handle] svg[y-slider] rect[pan-handle] circle[min-handle] circle[max-handle] canvas[canvas] canvas[overlay] i[icon-font-icon-settings] i[icon-plus] i[icon-minus] i[icon-font-icon-menu] div[viewer-Grid] svg[x-slider] rect[pan-handle] circle[min-handle] circle[max-handle] svg[y-slider] rect[pan-handle] circle[min-handle] circle[max-handle] canvas[canvas] canvas[overlay] i[icon-font-icon-settings] i[icon-plus] i[icon-minus] i[icon-font-icon-menu] rect[pan-handle] circle[min-handle] circle[max-handle]
dock tabs: 1 | Toolbox | Browse | Toolbox | Table | Table | Table | Table | MCL | Sequence Variability Map | Most Potent Residues | Logo Summary Table | Selection
model settings: {"settings":{"sequenceColumnName":"AlignedSequence","activityColumnName":"IC50","activityScaling":"none","columns":{},"showDendrogram":false,"showSequenceSpace":false,"sequenceSpaceParams":{"distanceF":"Needlemann-Wunsch","gapOpen":1.5,"gapExtend":0.2,"clusterEmbeddings":true,"epsilon":0.01,"minPts":4,"fingerprintType":"Morgan"},"mclSettings":{"maxIterations":16,"inflation":1.4,"threshold":70,"distanceF":"Needlemann-Wunsch","gapOpen":1.5,"gapExtend":0.2,"fingerprintType":"Morgan","useWebGPU":false,"minClusterSize":5,"webGPUDescription":"WebGPU is not supported on this device","webGPUDescriptionPromise":{}}},"cols":["Activity","AlignedSequence","ID","IC50","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","EmbedX (MCL)","EmbedY (MCL)","Cluster (MCL)","Cluster size (MCL)","Connectivity (MCL)"],"mcl":["EmbedX (MCL)","EmbedY (MCL)","Cluster (MCL)","Cluster size (MCL)","Connectivity (MCL)"],"seqSpace":[]}
svm grid status sample: {"areas":[["cell 2 of AAR",{"x":28,"y":20,"width":40,"height":20}],["cell 3 of AAR",{"x":28,"y":40,"width":40,"height":20}],["cell 22 of AAR",{"x":28,"y":60,"width":40,"height":20}],["cell 4 of AAR",{"x":28,"y":80,"width":40,"height":20}],["cell 5 of AAR",{"x":28,"y":100,"width":40,"height":20}]],"vals":[["rows shown",22],["rows",22],["pinned rows",0],["sort column","AAR"],["sort direction","ascending"],["is heatmap",false],["heatmap colors",true],["global color scaling",false],["max heatmap columns",100],["row height",20],["col labels orientation","Auto"],["effective col labels orientation","Horz"],["x scroll span",0.39861111111111114],["y scroll span",0.8522727272727273],["columns shown",9],["column order","AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"],["row order","2, 3, 22, 4, 5, 6, 20, 7, 8, 9"],["current row",2],["current column","∑"]],"textSample":[["text of cell 2 of AAR","A"],["text of cell 3 of AAR","C"],["text of cell 22 of AAR","COOH"],["text of cell 4 of AAR","D"],["text of cell 5 of AAR","E"],["text of cell 6 of AAR","F"],["text of cell 20 of AAR","G"],["text of cell 7 of AAR","H"]]}

```
