
function toParString(strings: TemplateStringsArray, values: any[]) {
  if (strings.length === 1)
    return strings.raw[0];
  else {
    let r = '';
    for (let i = 0; i < strings.length; i++)
      r += strings[i] + (values[i] ?? '');

    return r;
  }
}

const css = function(strings: TemplateStringsArray, ...values: any[]): CSSStyleSheet {
  const cssStyleSheet = new CSSStyleSheet();
  cssStyleSheet.replaceSync(toParString(strings, values));
  return cssStyleSheet;
};

export const style = css`
    :host {
        display: block;
    }`;

export const style1 = css`
/************* Panel with title bar ************/
.panel-base {
}
.panel-titlebar {
  flex-shrink: 0;
  background-color: #ffffff;
  color: #9497a0;
  height: 24px;
  width: 100%;
  overflow: hidden;
  padding: 0 2px;
  display: flex;
  align-items: center;
  box-sizing: border-box;
}
.panel-titlebar:hover {
  background-color: #f2f2f5;
  color: #9497a0;
}
.panel-titlebar-text {
  color: #9497a0;
}
.dockspan-tab-handle:only-child .dockspan-tab-handle-text {
  color: #9497a0;
}
.dockspan-tab-handle:only-child {
  border-right: none;
}
.dockspan-tab-handle:only-child .dockspan-tab-handle-text {
  margin-left: 2px;
  font-weight: normal;
}
.panel-titlebar-button-close {
  color:#aaa;
}
.panel-titlebar-button-close:hover {
  color:red;
}
.panel-content {
  /*background-color: #FFF;*/
    overflow: auto;
}
.panel-content * {
}
.panel-element-content-container {
  background-color: white;
  z-index: 1;
}
/***************** Floating dialog box ****************/
.dialog-floating {
  box-shadow: 5px 5px 20px #000;
  pointer-events: none;
}
/************ Dragging decorator ************/
.draggable-dragging-active {
  opacity: 0.95;
}
/************ Resize decorator ************/
.resize-handle {
  pointer-events: auto;
}
.resize-handle-corner {
  pointer-events: auto;
}
.resize-handle-e  {  }
.resize-handle-w  {  }
.resize-handle-s  {  }
.resize-handle-n  {  }
.resize-handle-ne {  }
.resize-handle-nw {  }
.resize-handle-se {  }
.resize-handle-sw {  }
/******************* Dock Manager ********************/
.dock-container {
  background-color: #888;
}
.dock-container-fill {
}
/******************* Document Manager ********************/
.document-manager {
  background-color: #666;
}
/**************************** Splitter *********************************/
.splitbar-horizontal {
  background-color: #f2f2f5;
  box-sizing: content-box;
  height: 1px;
  display: flex;
  padding: 0 2px;
  margin: 0 -2px;
  background-clip: content-box;
  width: 100% !important;
  float: left;
  z-index: 100;
  cursor: ns-resize;
  flex-shrink: 0;
  max-height: 1px;
}
.splitbar-horizontal:hover {
  background-color: #50A9C5;
}
.splitbar-vertical {
  background-color: #f2f2f5;
  box-sizing: content-box;
  width: 1px;
  display: flex;
  padding: 0 2px;
  margin: 0 -2px;
  background-clip: content-box;
  height: 100% !important;
  float: left;
  z-index: 100;
  cursor: ew-resize;
  flex-shrink: 0;
  max-width: 1px;
}
.splitbar-vertical:hover {
  background-color: #50A9C5;
}
.splitbar-horizontal-ghoust{
  background-color: #ffcc00;
}
.splitbar-vertical-ghoust {
  background-color: #ffcc00;
}
/*************************** Tab Host ********************************/
.dockspan-tab-content {
  background-color: #fff;
}
.dockspan-tab-content > * {
  background-color: #fff;
}
/*.tab-content * {
    margin: 0px;
}*/
.dockspan-tab-handle {
  background-color: #ffffff;
}
.dockspan-tab-handle:hover {
  background-color: #f2f2f5;
  color: #9497a0;
}
.dockspan-tab-handle-text {
}
.dockspan-tab-handle-close-button {
}
.dockspan-tab-handle-close-button:hover {
  color: red;
}
.dockspan-tab-handle-list-container {
  background-color: #ffffff;
}
.dockspan-tab-handle-content-seperator {
  background-color:  #333;
}
.dockspan-tab-handle-content-seperator-active {
  background-color: #008749;
}
`;

export const style2 = css`/************* Panel with title bar ************/
.panel-base {
  pointer-events: none;
}
.panel-base:focus {
  outline: 0;
}
.panel-titlebar, .dockspan-tab-handle:only-child{
  width: 100%;
  overflow: hidden;
  height: 25px;
  pointer-events: auto;
  z-index: 1;
  position: relative;
}
.panel-titlebar-icon {
  height: 16px;
  margin-right: 3px;
}
.panel-titlebar-button-close {
  visibility: hidden;
  margin-left: auto;
  margin-right: 3px;
  color: #9497a0;
  display: flex;
  justify-content: center;
  align-items: center;
}
.dockspan-tab-handle:only-child::after .dockspan-tab-handle-close-button {
  visibility: hidden;
  margin-left: auto;
  margin-right: 3px;
  color: #9497a0;
  display: flex;
  justify-content: center;
  align-items: center;
}
.panel-titlebar:hover .panel-titlebar-button-close {
  visibility: visible;
}
.dockspan-tab-handle:hover:only-child::after .dockspan-tab-handle-close-button {
   visibility: visible;
}
.panel-titlebar-button-close::after {
  float: right;
  cursor: pointer;
  content: "\e900";
  font-weight: normal;
  font-size: 18px;
  font-family: "icomoon";
  display: flex;
  height: 100%;
  align-items: center;
}
.dockspan-tab-handle:only-child .dockspan-tab-handle-close-button::after {
  float: right;
  cursor: pointer;
  content: "\e900";
  font-weight: normal;
  font-size: 18px;
  font-family: "icomoon";
  display: flex;
  height: 100%;
  align-items: center;
}
.panel-titlebar-button-close:hover {
  color: black;
}
.panel-content {
  width: 100%;
  overflow: hidden;
}
.panel-content * {
  /* margin: 0px; */
}
.panel-content-wrapper {
  width: 100%;
  height: 100%;
  position: relative;
}
.panel-grayout {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background: rgb(0 0 0 / 50%);
  z-index: 10;
  pointer-events: all;
}
.panel-has-changes {
    font-style: italic;
}
.panel-has-changes::before {
    content: '* ';
}
/***************** Floating dialog box ****************/
.dialog-floating {
  position: absolute;
}
.dialog-floating:focus {
  outline: 0;
}
/************ Resize decorator ************/
.resize-handle {
  position: absolute;
  width: 6px;
  height: 6px;
  z-index: 10;
}
.resize-handle-corner {
  position: absolute;
  width: 12px;
  height: 12px;
  z-index: 10;
}
.resize-handle-e {
  cursor: e-resize;
}
.resize-handle-w {
  cursor: w-resize;
}
.resize-handle-s {
  cursor: s-resize;
}
.resize-handle-n {
  cursor: n-resize;
}
.resize-handle-ne {
  cursor: ne-resize;
}
.resize-handle-nw {
  cursor: nw-resize;
}
.resize-handle-se {
  cursor: se-resize;
}
.resize-handle-sw {
  cursor: sw-resize;
}
/******************* Dock Manager ********************/
.dock-container {
  position: relative;
}
.dock-container-fill {}
/******************* Document Manager ********************/
.document-manager {}
/******************* Dock Wheel ********************/
.dock-wheel-base {
  position: absolute;
}
.dock-wheel-item {
  position: absolute;
  width: 32px;
  height: 32px;
}
.dock-wheel-fill {
  margin-left: -16px;
  margin-top: -16px;
}
.dock-wheel-left {
  margin-left: -48px;
  margin-top: -16px;
}
.dock-wheel-right {
  margin-left: 16px;
  margin-top: -16px;
}
.dock-wheel-top {
  margin-left: -16px;
  margin-top: -48px;
}
.dock-wheel-down {
  margin-left: -16px;
  margin-top: 16px;
}
.dock-wheel-panel-preview {
  position: absolute;
  background-color: rgba(128, 128, 255, 0.5);
}
.dock-wheel-fill-icon {
  background: url(../images/dock_fill.png) 0 0;
}
.dock-wheel-left-icon {
  background: url(../images/dock_left.png) 0 0;
}
.dock-wheel-right-icon {
  background: url(../images/dock_right.png) 0 0;
}
.dock-wheel-top-icon {
  background: url(../images/dock_top.png) 0 0;
}
.dock-wheel-down-icon {
  background: url(../images/dock_bottom.png) 0 0;
}
.dock-wheel-fill-icon-hover {
  background: url(../images/dock_fill_sel.png) 0 0;
}
.dock-wheel-left-icon-hover {
  background: url(../images/dock_left_sel.png) 0 0;
}
.dock-wheel-right-icon-hover {
  background: url(../images/dock_right_sel.png) 0 0;
}
.dock-wheel-top-icon-hover {
  background: url(../images/dock_top_sel.png) 0 0;
}
.dock-wheel-down-icon-hover {
  background: url(../images/dock_bottom_sel.png) 0 0;
}
/**************************** Splitter *********************************/
.splitter-container-horizontal {
  float: left;
  position: relative;
  pointer-events: none;
}
.splitter-container-vertical {
  position: relative;
  pointer-events: none;
}
.splitbar-horizontal {
  width: 100%;
  height: 5px;
  cursor: n-resize;
  position: relative;
  pointer-events: auto;
}
.splitbar-vertical {
  width: 5px;
  height: 100%;
  float: left;
  cursor: e-resize;
  position: relative;
  pointer-events: auto;
}
.splitbar-horizontal-ghoust {
  width: 100%;
  height: 5px;
  cursor: n-resize;
  position: absolute;
}
.splitbar-vertical-ghoust {
  width: 5px;
  height: 100%;
  position: absolute;
  cursor: e-resize;
}
/*************************** Tab Host ********************************/
.dockspan-tab-host {
  display: inline-block;
  position: absolute;
  background-color: #fff;
  left: 0px;
}
.dockspan-tab-content {
  position: relative;
}
.dockspan-tab-content:focus {
  outline: 0;
}
.dockspan-tab-content * {}
.dockspan-tab-handle {
  position: relative;
  height: 22px;
  float: left;
  overflow: hidden;
  cursor: pointer;
  padding-right: 16px;
  border-right: 1px solid var(--grey-2)
}
.dockspan-tab-handle:hover .dockspan-tab-handle-close-button {
  cursor: pointer;
  visibility: visible;
}
.dockspan-tab-handle-selected .dockspan-tab-handle-text {
  font-weight: bold;
  color: var(--grey-6)
}
.dockspan-tab-handle-text {
  margin-top: 3px;
  margin-left: 6px;
  margin-right: 6px;
  white-space: nowrap;
  overflow: hidden;
  float: left;
  text-overflow: ellipsis;
  width: calc(100% - 3px);
  display: block;
  color: #9497a0
}
.dockspan-tab-handle-text > span {
  vertical-align: top;
}
.dockspan-tab-handle-close-button {
  visibility: hidden;
}
.dockspan-tab-handle-close-button::after {
  float: right;
  position: absolute;
  right: 5px;
  content: "\e900";
  font-weight: normal;
  font-family: "icomoon";
  cursor: pointer;
  display: flex;
  height: 100%;
  align-items: center;
}
.dockspan-tab-handle-list-container {
  height: 22px;
  overflow: hidden;
  display: none;
  pointer-events: auto;
}
.dockspan-tab-handle-list-container-visible {
  display: flex;
}
.dockspan-tab-handle-content-seperator {
  height: 4px;
  display: none;
}
.dockspan-tab-handle-content-seperator-visible {
  display: none;
}
.dockspan-tab-handle-content-seperator-selected {}
.dockspab-tab-handle-context-menu {
  position: absolute;
  z-index: 10000000;
  background-color: white;
  border: white;
  border: black solid 2px;
  padding: 5px;
}
.dockspab-tab-handle-context-menu>div {
  cursor: pointer;
  padding: 2px;
}
.dockspab-tab-handle-context-menu>div:hover {
  background-color: gray;
  cursor: pointer;
  padding: 2px;
}
/*************************** Text Selection **************************/
.disable-selection {
  user-select: none;
  -moz-user-select: none;
  -webkit-user-select: none;
  -ms-user-select: none;
  cursor: default;
}`;
