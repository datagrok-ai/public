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
import {FlowNode, FlowScheme, EXEC_IN_KEY, EXEC_OUT_KEY, ORDER_SOCKET_TYPE, isExecKey, nodeMissingRequirements} from './scheme';
import {TypedSocket} from './sockets';
import {getSlotColor, getSlotLetter, pastelize} from '../types/type-map';
import {tid} from '../utils/test-ids';
import {summarizeNode} from '../summary/summary-generator';

interface NodeProps {
  data: FlowNode & {selected?: boolean};
  emit: RenderEmit<FlowScheme>;
}

/** The host `FlowEditor` exposes a narrow callback surface on every node it
 *  owns (`FlowNode.editorBridge`, stamped when the node enters the editor's
 *  data layer). Resolving it from the node — never from a global — keeps each
 *  component bound to its own editor: several editors coexist on a page (file
 *  previews, the creation-script dialog, detached compile editors). */
function toggleCollapsed(node: FlowNode): void {
  node.editorBridge?.toggleCollapsed(node.id);
}
function isConnected(node: FlowNode, side: 'input' | 'output', key: string): boolean {
  return node.editorBridge?.isSocketConnected(node.id, side, key) ?? false;
}

export function FlowNodeComponent(props: NodeProps): React.JSX.Element {
  const node = props.data;
  // Output nodes live docked in the Outputs strip (FlowEditor pins their
  // position) and render as compact rows, not full node cards.
  if (node.dgNodeType === 'output') return OutputRowComponent(props);
  const collapsed = node.collapsed === true;
  // Exec (execution-ordering) ports render separately at the top corners — keep
  // them out of the regular data-socket rows. Hidden inputs (and their
  // pass-throughs) render only when connected — their slots stay data-carrying,
  // the user just never sees an unwired one.
  const hiddenRow = (key: string, side: 'input' | 'output'): boolean => {
    const base = side === 'input' ? key : (key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : '');
    return node.hiddenInputs.has(base) && !isConnected(node, side, key);
  };
  const inputs = (Object.entries(node.inputs) as Array<[string, ClassicPreset.Input<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key) && !hiddenRow(key, 'input'));
  const outputs = (Object.entries(node.outputs) as Array<[string, ClassicPreset.Output<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key) && !hiddenRow(key, 'output'));
  const execIn = node.inputs[EXEC_IN_KEY] as ClassicPreset.Input<TypedSocket> | undefined;
  const execOut = node.outputs[EXEC_OUT_KEY] as ClassicPreset.Output<TypedSocket> | undefined;
  const controls = Object.entries(node.controls) as Array<[string, ClassicPreset.Control | undefined]>;
  const ptCount = node.passthroughCount ?? 0;

  // Title bars get the pastel of the node's identity color — vivid hues are
  // kept for small surfaces (minimap, sockets) where saturation aids legibility.
  const titleColor = (node as unknown as {color?: string}).color;
  const titleStyle: React.CSSProperties = titleColor ? {background: pastelize(titleColor)} : {};

  const dgStatus = (node as unknown as {dgStatus?: string}).dgStatus ?? 'idle';
  const statusText = (node as unknown as {statusText?: string}).statusText ?? '';

  // Pre-run hint: structural inputs the user still has to provide. Shown only
  // when the node hasn't successfully run (idle/stale), so a "Done"/"Error"
  // status from a real run always takes precedence.
  const needs = nodeMissingRequirements(node, (key) => isConnected(node, 'input', key));
  const idle = !dgStatus || dgStatus === 'idle' || dgStatus === 'stale';
  const attention = idle && needs.length > 0;

  // Auto-summary caption (U12), shown when the user hasn't written their own
  // description. Computed once so the visible (CSS-ellipsized) text and the
  // hover tooltip stay in sync — the tooltip reveals it in full when truncated.
  const autoSummary = !node.description && !collapsed ? summarizeNode(node) : '';

  // Collapse is the caret's job now — the status dot is display-only, so
  // clicking the run indicator never hides the node out from under you.
  const onCaretClick = (e: React.MouseEvent): void => {
    e.stopPropagation();
    toggleCollapsed(node);
  };
  const stopPointer = (e: React.PointerEvent): void => {
    // Keep the AreaPlugin's node-drag handler from grabbing a click on a
    // title-bar affordance (caret).
    e.stopPropagation();
  };

  return (
    <div
      className={`ff-node ff-node-${node.dgNodeType ?? 'func'}` + (collapsed ? ' ff-node-collapsed' : '')}
      data-testid={tid('node')}
      data-node-id={node.id}
      data-node-type={node.dgNodeType ?? 'func'}
      data-node-type-name={node.dgTypeName ?? ''}
      data-func={node.dgFuncName ?? ''}
      data-node-label={node.label}
      data-selected={node.selected ? 'true' : 'false'}
      data-status={dgStatus}
      data-attention={attention ? 'true' : 'false'}
    >
      {/* Execution-ordering ports — top corners (KNIME flow-variable style).
          exec-in (left) accepts "run after" predecessors; exec-out (right)
          drives successors. Always rendered so edges stay attached even when
          the node is collapsed — but only *visible* when one of them is wired,
          the node is hovered, or an order drag is in progress (CSS keys off
          data-wired / .ff-node:hover / .ff-connecting-order). */}
      <div
        className="ff-node-exec-row"
        data-wired={(isConnected(node, 'input', EXEC_IN_KEY) || isConnected(node, 'output', EXEC_OUT_KEY)) ?
          'true' : 'false'}
      >
        <span
          className="ff-exec-port ff-exec-in" data-testid={tid('exec-in')}
          title={'Run order: this node waits for its predecessors. Drag a wire from another node\'s ' +
            'order square to here to make this node run after it — sequencing only, no data flows.'}
        >
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
        <span
          className="ff-exec-port ff-exec-out" data-testid={tid('exec-out')}
          title={'Run order: drag a wire from here to another node\'s order square to make that node ' +
            'run after this one — sequencing only, no data flows.'}
        >
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

      <div className="ff-node-title" data-testid={tid('node-title')} style={titleStyle} data-role={node.dgRole ?? ''}>
        <div
          className="ff-node-status"
          data-testid={tid('node-status')}
          data-status={dgStatus}
          title={statusText || 'Not run yet'}
        />
        <span className="ff-node-title-text" data-testid={tid('node-title-text')}>{node.label}</span>
        <span
          className="ff-node-caret"
          data-testid={tid('node-caret')}
          title={collapsed ? 'Expand node' : 'Collapse node'}
          onPointerDown={stopPointer}
          onClick={onCaretClick}
        >{collapsed ? '▸' : '▾'}</span>
      </div>

      {/* Pre-run "Needs input" hint takes precedence over an idle/stale status;
          a real run's status (Done/Running/Error) wins otherwise. Shown
          collapsed too, so a folded node still reports its state at a glance. */}
      {attention ? (
        <div className="ff-node-hint" data-testid={tid('node-hint')} title={`Connect or set: ${needs.join(', ')}`}>
          Requires: {needs.join(', ')}
        </div>
      ) : statusText ? (
        <div className="ff-node-statusline" data-testid={tid('node-statusline')} data-status={dgStatus}>{statusText}</div>
      ) : null}

      {node.description && !collapsed && (
        <div className="ff-node-description" data-testid={tid('node-description')} title={node.description}>{node.description}</div>
      )}

      {/* Auto-generated plain-language caption — the flow documents itself
          (U12). The title reveals the full text when CSS ellipsis truncates it. */}
      {autoSummary && (
        <div className="ff-node-summary" data-testid={tid('node-summary')} title={autoSummary}>{autoSummary}</div>
      )}

      {!collapsed && (
        <div className="ff-node-body" data-testid={tid('node-body')}>
          <div className="ff-node-io">
            <div className="ff-node-inputs">
              {inputs.map(([key, input]) => input && (
                <div key={key} className="ff-socket-row ff-socket-row-input" data-testid={tid('socket-input', key)}
                  title={node.inputDescriptions?.[key] || undefined}>
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
              {outputs.map(([key, output]) => output && (
                <div
                  key={key}
                  className={
                    'ff-socket-row ff-socket-row-output' +
                    (ptCount > 0 && key.endsWith('__pt') ? ' ff-socket-row-passthrough' : '')
                  }
                  data-testid={tid('socket-output', key)}
                  title={node.outputDescriptions?.[key] || undefined}
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
            {inputs.map(([key, input]) => input && isConnected(node, 'input', key) && (
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
            {outputs.map(([key, output]) => output && isConnected(node, 'output', key) && (
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

/** Output nodes have no canvas presence at all — their visible form is a
 *  screen-space chip inside the Outputs strip (built by `FlowEditor`), and
 *  their wire endpoints come from the analytic chip-aware socket-position
 *  watcher, not from DOM. Render an empty hidden placeholder so the rete node
 *  view exists but never paints, scales, or intercepts pointer events. */
function OutputRowComponent(props: NodeProps): React.JSX.Element {
  return <div style={{display: 'none'}} data-node-id={props.data.id} />;
}

interface SocketProps {
  data: TypedSocket;
}

/** Socket chip: a single type letter colored like the Column Manager's type
 *  column (`t` table, `s` string, `i` int, `d` double, …; see `getSlotLetter`).
 *  The type color travels as the `--socket-color` CSS var — `.ff-socket` paints
 *  the letter and ring with it, and `.ff-socket-row-passthrough .ff-socket`
 *  reads it for the dashed ring. An `order` socket stays a bare gray square
 *  (no letter) and gets NO title of its own — a nested title would shadow the
 *  exec-port wrapper's plain-language "Run order: drag…" explanation with a
 *  bare "order". */
export function FlowSocketComponent(props: SocketProps): React.JSX.Element {
  const color = getSlotColor(props.data.dgType);
  const isOrder = props.data.dgType === ORDER_SOCKET_TYPE;
  return (
    <div
      className="ff-socket"
      data-testid={tid('socket', props.data.dgType)}
      style={{['--socket-color' as never]: color}}
      title={isOrder ? undefined : props.data.dgType}
      data-type={props.data.dgType}
    >{isOrder ? null : getSlotLetter(props.data.dgType)}</div>
  );
}

/** Custom Connection component that paints the path in the source slot's
 *  color. Replaces the default `Connection` that hardcodes `stroke: steelblue`
 *  via `styled-components` (whose CSS otherwise wins over a plain
 *  `setAttribute('stroke', ...)`). The color is stuffed into the connection
 *  payload by `FlowEditor` at construction time as `_color`. */
interface ConnectionProps {
  data: FlowScheme['Connection'] & {_color?: string; _count?: string};
}

export function FlowConnectionComponent(props: ConnectionProps): React.JSX.Element | null {
  const {path, start, end} = Presets.classic.useConnection();
  if (!path) return null;
  const color = props.data._color ?? '#8892a0';
  const count = props.data._count;

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
      data-testid={tid('connection')}
      data-connection-id={props.data.id}
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
          data-testid={tid('waypoint', i)}
          data-connection-id={props.data.id}
          data-waypoint-index={i}
          style={{fill: color, stroke: '#fff', strokeWidth: 1.5, pointerEvents: 'auto', cursor: 'move'}}
        />
      ))}
      {count && start && end && (
        <text
          className="ff-edge-count"
          data-testid={tid('edge-count')}
          x={(start.x + end.x) / 2}
          y={(start.y + end.y) / 2 - 6}
          textAnchor="middle"
          dominantBaseline="central"
        >{count}</text>
      )}
    </svg>
  );
}
