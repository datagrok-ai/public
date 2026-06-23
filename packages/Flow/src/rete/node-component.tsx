/** Custom React node component used by `rete-react-plugin`'s classic preset.
 *
 * Despite being React, the output is just plain DOM — no React state outside
 * this file, no React in the rest of the package. The component is mounted
 * by ReactPlugin into the empty `<div>` AreaPlugin creates per node, exactly
 * the same way `KetcherSketcher` mounts the Ketcher editor.
 *
 * Visual contract: the DOM structure here matches the `ff-*` selectors in
 * `css/funcflow.css`. Don't change classnames without updating the CSS. */

import * as React from 'react';
import {ClassicPreset} from 'rete';
import {Presets} from 'rete-react-plugin';
import type {RenderEmit} from 'rete-react-plugin';
import {classicConnectionPath} from 'rete-render-utils';

const {RefSocket, RefControl} = Presets.classic;
import {FlowNode, FlowScheme, EXEC_IN_KEY, EXEC_OUT_KEY, isExecKey} from './scheme';
import {TypedSocket} from './sockets';
import {getSlotColor} from '../types/type-map';

interface NodeProps {
  data: FlowNode & {selected?: boolean};
  emit: RenderEmit<FlowScheme>;
}

/** The host `FlowEditor` exposes a narrow callback surface on
 *  `window.__ff_editor` so React components can call back into it without
 *  threading the editor through React context (we have only one editor at a
 *  time per page). */
interface EditorBridge {
  toggleCollapsed(id: string): void;
  isSocketConnected(nodeId: string, side: 'input' | 'output', key: string): boolean;
}
function bridge(): EditorBridge | undefined {
  return (window as unknown as {__ff_editor?: EditorBridge}).__ff_editor;
}
function toggleCollapsed(nodeId: string): void {
  bridge()?.toggleCollapsed(nodeId);
}
function isConnected(nodeId: string, side: 'input' | 'output', key: string): boolean {
  return bridge()?.isSocketConnected(nodeId, side, key) ?? false;
}

export function FlowNodeComponent(props: NodeProps): React.JSX.Element {
  const node = props.data;
  const collapsed = node.collapsed === true;
  // Exec (execution-ordering) ports render separately at the top corners — keep
  // them out of the regular data-socket rows.
  const inputs = (Object.entries(node.inputs) as Array<[string, ClassicPreset.Input<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key));
  const outputs = (Object.entries(node.outputs) as Array<[string, ClassicPreset.Output<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key));
  const execIn = node.inputs[EXEC_IN_KEY] as ClassicPreset.Input<TypedSocket> | undefined;
  const execOut = node.outputs[EXEC_OUT_KEY] as ClassicPreset.Output<TypedSocket> | undefined;
  const controls = Object.entries(node.controls) as Array<[string, ClassicPreset.Control | undefined]>;
  const ptCount = node.passthroughCount ?? 0;

  const titleColor = (node as unknown as {color?: string}).color;
  const titleStyle: React.CSSProperties = titleColor ? {background: titleColor} : {};

  const onStatusClick = (e: React.MouseEvent): void => {
    e.stopPropagation();
    toggleCollapsed(node.id);
  };
  const onStatusPointerDown = (e: React.PointerEvent): void => {
    // Prevent the AreaPlugin's drag handler from picking up the pointerdown
    // when the user is just clicking the status dot to collapse.
    e.stopPropagation();
  };

  return (
    <div
      className={`ff-node ff-node-${node.dgNodeType ?? 'func'}` + (collapsed ? ' ff-node-collapsed' : '')}
      data-node-id={node.id}
      data-selected={node.selected ? 'true' : 'false'}
      data-status={(node as unknown as {dgStatus?: string}).dgStatus ?? 'idle'}
    >
      {/* Execution-ordering ports — top corners (KNIME flow-variable style).
          exec-in (left) accepts "run after" predecessors; exec-out (right)
          drives successors. Always rendered so edges stay attached even when
          the node is collapsed. */}
      <div className="ff-node-exec-row">
        <span className="ff-exec-port ff-exec-in" title="Run after (order in)">
          {execIn && (
            <RefSocket
              name="exec-in-socket"
              emit={props.emit}
              side="input"
              socketKey={EXEC_IN_KEY}
              nodeId={node.id}
              payload={execIn.socket}
            />
          )}
        </span>
        <span className="ff-exec-port ff-exec-out" title="Run before (order out)">
          {execOut && (
            <RefSocket
              name="exec-out-socket"
              emit={props.emit}
              side="output"
              socketKey={EXEC_OUT_KEY}
              nodeId={node.id}
              payload={execOut.socket}
            />
          )}
        </span>
      </div>

      <div className="ff-node-title" style={titleStyle} data-role={node.dgRole ?? ''}>
        <div
          className="ff-node-status"
          data-status={(node as unknown as {dgStatus?: string}).dgStatus ?? 'idle'}
          title={collapsed ? 'Expand node' : 'Collapse node'}
          onPointerDown={onStatusPointerDown}
          onClick={onStatusClick}
        />
        <span className="ff-node-title-text">{node.label}</span>
      </div>

      {node.description && !collapsed && (
        <div className="ff-node-description" title={node.description}>{node.description}</div>
      )}

      {!collapsed && (
        <div className="ff-node-body">
          <div className="ff-node-io">
            <div className="ff-node-inputs">
              {inputs.map(([key, input]) => input && (
                <div key={key} className="ff-socket-row ff-socket-row-input">
                  <RefSocket
                    name="input-socket"
                    emit={props.emit}
                    side="input"
                    socketKey={key}
                    nodeId={node.id}
                    payload={input.socket}
                  />
                  <span className="ff-socket-label">{input.label ?? key}</span>
                </div>
              ))}
            </div>

            <div className="ff-node-outputs">
              {outputs.map(([key, output], idx) => output && (
                <div
                  key={key}
                  className={
                    'ff-socket-row ff-socket-row-output' +
                    (idx < ptCount ? ' ff-socket-row-passthrough' : '')
                  }
                >
                  <span className="ff-socket-label">{output.label ?? key}</span>
                  <RefSocket
                    name="output-socket"
                    emit={props.emit}
                    side="output"
                    socketKey={key}
                    nodeId={node.id}
                    payload={output.socket}
                  />
                </div>
              ))}
            </div>
          </div>

          {controls.length > 0 && (
            <div className="ff-node-controls">
              {controls.map(([key, control]) => control && (
                <RefControl key={key} name="control" emit={props.emit} payload={control} />
              ))}
            </div>
          )}
        </div>
      )}

      {/* Collapsed nodes only render sockets that have at least one connection.
          Disconnected sockets are useless without the body (you can't drag a
          new wire from a tiny dot you can barely see), so showing them is just
          clutter. Connected ones still need DOM so existing wires keep their
          endpoints. */}
      {collapsed && (
        <div className="ff-node-collapsed-sockets">
          <div className="ff-collapsed-inputs">
            {inputs.map(([key, input]) => input && isConnected(node.id, 'input', key) && (
              <RefSocket
                key={key}
                name="input-socket"
                emit={props.emit}
                side="input"
                socketKey={key}
                nodeId={node.id}
                payload={input.socket}
              />
            ))}
          </div>
          <div className="ff-collapsed-outputs">
            {outputs.map(([key, output]) => output && isConnected(node.id, 'output', key) && (
              <RefSocket
                key={key}
                name="output-socket"
                emit={props.emit}
                side="output"
                socketKey={key}
                nodeId={node.id}
                payload={output.socket}
              />
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

interface SocketProps {
  data: TypedSocket;
}

/** Colored socket dot. We expose the type color twice:
 *  - `background` paints the regular (non-pass-through) dot
 *  - `--socket-color` CSS var is read by `.ff-socket-row-passthrough .ff-socket`
 *    to color the dashed ring (where `background` is forced white). */
export function FlowSocketComponent(props: SocketProps): React.JSX.Element {
  const color = getSlotColor(props.data.dgType);
  return (
    <div
      className="ff-socket"
      style={{background: color, ['--socket-color' as never]: color}}
      title={props.data.dgType}
      data-type={props.data.dgType}
    />
  );
}

/** Custom Connection component that paints the path in the source slot's
 *  color. Replaces the default `Connection` that hardcodes `stroke: steelblue`
 *  via `styled-components` (whose CSS otherwise wins over a plain
 *  `setAttribute('stroke', ...)`). The color is stuffed into the connection
 *  payload by `FlowEditor` at construction time as `_color`. */
interface ConnectionProps {
  data: FlowScheme['Connection'] & {_color?: string};
}

export function FlowConnectionComponent(props: ConnectionProps): React.JSX.Element | null {
  const {path, start, end} = Presets.classic.useConnection();
  if (!path) return null;
  const color = props.data._color ?? '#8892a0';

  // Compute the actual path: if waypoints are present, chain
  // classicConnectionPath segments through start → waypoints → end. Without
  // waypoints, fall back to the preset's path (computed by useConnection).
  const waypoints = (props.data as {waypoints?: Array<{x: number; y: number}>}).waypoints;
  let drawPath = path;
  if (waypoints && waypoints.length > 0 && start && end) {
    const points = [start, ...waypoints, end];
    drawPath = '';
    for (let i = 0; i < points.length - 1; i++)
      drawPath += classicConnectionPath([points[i], points[i + 1]], 0.3);
  }

  // Two paths: a wide invisible "hit" path (transparent stroke, 14px) so
  // right-click and pointer events land reliably even though the visible
  // line is only 2.5px, and the visible path itself with pointerEvents:none
  // (the hit path handles all interaction). pointerEvents: 'stroke' on the
  // hit path means events fire only within the stroke area, not in the
  // unbounded SVG canvas.
  return (
    <svg
      data-testid="connection"
      style={{
        overflow: 'visible',
        position: 'absolute',
        pointerEvents: 'none',
        width: '9999px',
        height: '9999px',
      }}
    >
      <path
        d={drawPath}
        fill="none"
        stroke="transparent"
        strokeWidth={14}
        style={{pointerEvents: 'stroke'}}
      />
      <path
        d={drawPath}
        className="ff-connection-path"
        style={{stroke: color, pointerEvents: 'none'}}
      />
      {waypoints?.map((wp, i) => (
        <circle
          key={`wp-${i}`}
          cx={wp.x}
          cy={wp.y}
          r={5}
          className="ff-waypoint"
          data-connection-id={props.data.id}
          data-waypoint-index={i}
          style={{fill: color, stroke: '#fff', strokeWidth: 1.5, pointerEvents: 'auto', cursor: 'move'}}
        />
      ))}
    </svg>
  );
}
