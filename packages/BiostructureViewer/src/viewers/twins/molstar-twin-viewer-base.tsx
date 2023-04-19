import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as React from 'react';
import * as ReactDom from 'react-dom';
import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {CSSProperties} from 'react';
import {RcsbFv3DComponent, RcsbFv3DCssConfig} from '../../libs/rcsb-saguaro-3d/RcsbFv3D/RcsbFv3DComponent';
import {RcsbFvSequenceInterface} from '../../libs/rcsb-saguaro-3d/RcsbFvSequence/RcsbFvSequence';
import {EventType, RcsbFvContextManager} from '../../libs/rcsb-saguaro-3d/RcsbFvContextManager/RcsbFvContextManager';
import {RcsbFvStructureInterface} from '../../libs/rcsb-saguaro-3d/RcsbFvStructure/RcsbFvStructure';

export interface RcsbFv3DAbstractInterface {
  elementId: string;
  cssConfig?: RcsbFv3DCssConfig;
}

export abstract class RcsbFv3dBase {
  protected elementId: string;
  protected structureConfig: RcsbFvStructureInterface;
  protected sequenceConfig: RcsbFvSequenceInterface;
  protected ctxManager: RcsbFvContextManager = new RcsbFvContextManager();
  private fullScreenFlag: boolean = false;
  protected cssConfig: {
    rootPanel?: CSSProperties,
    structurePanel?: CSSProperties,
    sequencePanel?: CSSProperties
  } | undefined;

  constructor(config?: any) {
    if (config != null)
      this.init(config);
  }

  protected abstract init(config: any): void;

  public render(): void {
    if (this.elementId == null)
      throw new Error('HTML element not found');
    const element: HTMLElement = document.getElementById(this.elementId) ?? document.createElement<'div'>('div');
    if (element.getAttribute('id') == null) {
      element.setAttribute('id', this.elementId);
      document.body.append(element);
      this.fullScreenFlag = true;
      document.body.style.overflow = 'hidden';
    }
    ReactDom.render(
      <RcsbFv3DComponent
        structurePanelConfig={this.structureConfig}
        sequencePanelConfig={this.sequenceConfig}
        id={'RcsbFv3D_innerDiv_' + Math.random().toString(36).substr(2)}
        ctxManager={this.ctxManager}
        cssConfig={this.cssConfig}
        unmount={this.unmount.bind(this)}
        fullScreen={this.fullScreenFlag}
      />,
      element
    );
  }

  public unmount(removeHtmlElement?: boolean): void {
    const element: HTMLElement | null = document.getElementById(this.elementId);
    if (element != null) {
      ReactDom.unmountComponentAtNode(element);
      if (removeHtmlElement) {
        element.remove();
        document.body.style.overflow = 'visible';
      }
      window.history.back();
    }
  }

  public updateConfig(config: {
    structurePanelConfig?: RcsbFvStructureInterface; sequencePanelConfig?: RcsbFvSequenceInterface;
  }) {
    this.ctxManager.next({eventType: EventType.UPDATE_CONFIG, eventData: config});
  }

  public pluginCall(f: (plugin: PluginContext) => void) {
    this.ctxManager.next({eventType: EventType.PLUGIN_CALL, eventData: f});
  }
}
