/** React node components for rete-react-plugin. The DOM structure matches the
 *  `ff-*` selectors in `css/funcflow.css`. */

import * as React from 'react';
import * as ui from 'datagrok-api/ui';
import {ClassicPreset} from 'rete';
import {Presets} from 'rete-react-plugin';
import type {RenderEmit} from 'rete-react-plugin';
import {classicConnectionPath} from 'rete-render-utils';

const {RefSocket, RefControl} = Presets.classic;
import {FlowNode, FlowScheme, EXEC_IN_KEY, EXEC_OUT_KEY, ORDER_SOCKET_TYPE, isExecKey, nodeMissingRequirements,
  hiddenSocketRow, supportsInlinePreview, inlinePreviewEnabled, inlinePreviewSize,
  INLINE_PREVIEW_SIZE_PROP} from './scheme';
import {InputValueControl} from './nodes/input-value-control';
import {TypedSocket} from './sockets';
import {getSlotColor, getSlotLetter, pastelize} from '../types/type-map';
import {tid} from '../utils/test-ids';
import {summarizeNode} from '../summary/summary-generator';

interface NodeProps {
  data: FlowNode & {selected?: boolean};
  emit: RenderEmit<FlowScheme>;
}

// Resolve the bridge from the node, never a global — several editors coexist on a page.
function toggleCollapsed(node: FlowNode): void {
  node.editorBridge?.toggleCollapsed(node.id);
}
function isConnected(node: FlowNode, side: 'input' | 'output', key: string): boolean {
  return node.editorBridge?.isSocketConnected(node.id, side, key) ?? false;
}

export function FlowNodeComponent(props: NodeProps): React.JSX.Element {
  const node = props.data;
  if (node.dgNodeType === 'output') return OutputRowComponent(props);
  const collapsed = node.collapsed === true;
  const hiddenRow = (key: string, side: 'input' | 'output'): boolean =>
    hiddenSocketRow(node, side, key, (s, k) => isConnected(node, s, k));
  const inputs = (Object.entries(node.inputs) as Array<[string, ClassicPreset.Input<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key) && !hiddenRow(key, 'input'));
  const outputs = (Object.entries(node.outputs) as Array<[string, ClassicPreset.Output<TypedSocket> | undefined]>)
    .filter(([key]) => !isExecKey(key) && !hiddenRow(key, 'output'));
  const execIn = node.inputs[EXEC_IN_KEY] as ClassicPreset.Input<TypedSocket> | undefined;
  const execOut = node.outputs[EXEC_OUT_KEY] as ClassicPreset.Output<TypedSocket> | undefined;
  const controls = Object.entries(node.controls) as Array<[string, ClassicPreset.Control | undefined]>;
  const ptCount = node.passthroughCount ?? 0;

  const titleColor = (node as unknown as {color?: string}).color;
  const titleStyle: React.CSSProperties = titleColor ? {background: pastelize(titleColor)} : {};

  const dgStatus = (node as unknown as {dgStatus?: string}).dgStatus ?? 'idle';
  const statusText = (node as unknown as {statusText?: string}).statusText ?? '';

  const needs = nodeMissingRequirements(node, (key) => isConnected(node, 'input', key));
  const idle = !dgStatus || dgStatus === 'idle' || dgStatus === 'stale';
  const attention = idle && needs.length > 0;

  const autoSummary = !node.description && !collapsed ? summarizeNode(node) : '';

  const previewCapable = supportsInlinePreview(node);
  const previewOn = previewCapable && inlinePreviewEnabled(node);

  const onCaretClick = (e: React.MouseEvent): void => {
    e.stopPropagation();
    toggleCollapsed(node);
  };
  const stopPointer = (e: React.PointerEvent): void => {
    // keep AreaPlugin's node-drag handler off title-bar affordances
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
      data-inline-preview={previewOn ? 'true' : 'false'}
    >
      {/* Exec ports always render so edges keep their endpoints; CSS shows them
          only when wired, hovered, or during an order drag. */}
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
        <span className="ff-node-title-text" data-testid={tid('node-title-text')} title={node.label}>{node.label}</span>
        {previewCapable && (
          <span
            className="ff-node-preview-toggle"
            data-testid={tid('node-preview-toggle')}
            data-on={previewOn ? 'true' : 'false'}
            title={previewOn ? 'Hide the in-node preview' : 'Show the result right on the node'}
            onPointerDown={stopPointer}
            onClick={(e) => {
              e.stopPropagation();
              node.editorBridge?.toggleInlinePreview(node.id);
            }}
          >{previewOn ? '⊟' : '⊞'}</span>
        )}
        <span
          className="ff-node-caret"
          data-testid={tid('node-caret')}
          title={collapsed ? 'Expand node' : 'Collapse node'}
          onPointerDown={stopPointer}
          onClick={onCaretClick}
        >{collapsed ? '▸' : '▾'}</span>
      </div>

      {/* One always-rendered info line (hint > status > description > summary) —
          the card's height never changes across run states. */}
      {!collapsed && (() => {
        const line = attention ?
          {kind: 'hint', text: `Requires: ${needs.join(', ')}`, tip: `Connect or set: ${needs.join(', ')}`,
            testid: tid('node-hint')} :
          statusText ?
            {kind: 'status', text: statusText, tip: statusText, testid: tid('node-statusline')} :
            node.description ?
              {kind: 'description', text: node.description, tip: node.description, testid: tid('node-description')} :
              autoSummary ?
                {kind: 'summary', text: autoSummary, tip: autoSummary, testid: tid('node-summary')} :
                {kind: 'empty', text: ' ', tip: '', testid: tid('node-infoline')};
        return (
          <div
            className="ff-node-infoline" data-testid={line.testid} data-kind={line.kind}
            data-status={dgStatus} title={line.tip}
          >{line.text}</div>
        );
      })()}

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
              {/* ⋯ indicator: hidden toggleable inputs — click pops the "Shown inputs" checkboxes. */}
              {(() => {
                const hiddenCount = Object.keys(node.inputs).filter((k) =>
                  !isExecKey(k) && !node.hiddenInputs.has(k) && hiddenRow(k, 'input')).length;
                return hiddenCount > 0 && (
                  <div
                    className="ff-node-more-inputs"
                    data-testid={tid('node-more-inputs')}
                    ref={(el) => {
                      if (el) ui.tooltip.bind(el, `${hiddenCount} input${hiddenCount > 1 ? 's are' : ' is'} ` +
                        'hidden — click to choose which are shown');
                    }}
                    onPointerDown={stopPointer}
                    onClick={(e) => {
                      e.stopPropagation();
                      node.editorBridge?.showShownInputsMenu(node.id, e.nativeEvent);
                    }}
                  >⋯</div>
                );
              })()}
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

          {previewOn && <InlineNodePreview node={node} />}
        </div>
      )}

      {/* Collapsed nodes render socket DOM only for connected sockets, so wires keep endpoints. */}
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

/** In-node preview of a viewer/widget output. The in-card container only
 *  reserves the box (border, placeholder, resize) — the actual content is
 *  mounted by the editor in a SCREEN-SPACE PORTAL tracking this container
 *  (`FlowEditor.syncInlinePreview`): DG viewer popups are `position: fixed`,
 *  and any transformed ancestor (the zoom/pan canvas) would become their
 *  containing block, throwing them far from the viewer. The size persists in
 *  `node.properties` so it serializes with the flow. */
function InlineNodePreview(props: {node: FlowNode}): React.JSX.Element {
  const node = props.node;
  const hostRef = React.useRef<HTMLDivElement>(null);

  // Placeholder before content exists; the portal binds/refreshes on every
  // render, so a fresh run's content appears without extra plumbing.
  React.useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    const has = (node.editorBridge?.getInlinePreviewContent(node.id) ?? null) != null;
    if (has) {
      host.querySelector(':scope > .ff-node-preview-placeholder')?.remove();
      host.dataset.empty = 'false';
    } else {
      // A run on its way to this node shows a loader, not the resting hint —
      // a fresh full run clears the captured value before recomputing it.
      const pending = node.editorBridge?.isInlinePreviewPending(node.id) ?? false;
      const want = pending ? 'loading' : 'true';
      if (host.dataset.empty !== want) {
        host.innerHTML = '';
        const ph = document.createElement('div');
        ph.className = 'ff-node-preview-placeholder';
        ph.dataset.testid = tid('node-preview-placeholder');
        if (pending) {
          ph.dataset.loading = 'true';
          ph.appendChild(ui.loader());
        } else
          ph.textContent = 'Run the flow to see the preview';
        host.appendChild(ph);
        host.dataset.empty = want;
      }
    }
    node.editorBridge?.syncInlinePreview(node.id, host);
  });

  // Tear the portal down when the preview unmounts (toggle off, collapse, node
  // removed) — this also releases the hosted marker to the bottom panel.
  React.useEffect(() => () => {
    node.editorBridge?.releaseInlinePreview(node.id);
  }, []);

  // Persist a user drag-resize into the node properties (cosmetic — no
  // params-changed report; properties serialize on save).
  React.useEffect(() => {
    const host = hostRef.current;
    if (!host || typeof ResizeObserver === 'undefined') return;
    const ro = new ResizeObserver(() => {
      // The CSS `resize` handle writes inline width/height — read those, not
      // offsetWidth (which adds the borders and would drift the stored size).
      const w = Math.round(parseFloat(host.style.width) || host.offsetWidth);
      const h = Math.round(parseFloat(host.style.height) || host.offsetHeight);
      if (w <= 0 || h <= 0) return;
      const cur = inlinePreviewSize(node);
      if (cur.width !== w || cur.height !== h)
        node.properties[INLINE_PREVIEW_SIZE_PROP] = {width: w, height: h};
    });
    ro.observe(host);
    return () => ro.disconnect();
  }, []);

  const size = inlinePreviewSize(node);
  return (
    <div
      className="ff-node-inline-preview"
      data-testid={tid('node-preview')}
      style={{width: `${size.width}px`, height: `${size.height}px`}}
      ref={hostRef}
      onPointerDown={(e) => e.stopPropagation()}
      onDoubleClick={(e) => e.stopPropagation()}
      onWheel={(e) => e.stopPropagation()}
    />
  );
}

/** Hosts a DG input in a node body. The element is built once and re-attached
 *  on every (re)mount, so React re-renders never rebuild it or drop focus. */
export function DgControlComponent(props: {data: InputValueControl}): React.JSX.Element {
  const ref = React.useRef<HTMLDivElement>(null);
  React.useEffect(() => {
    const el = props.data.element();
    if (el && ref.current && el.parentElement !== ref.current)
      ref.current.appendChild(el);
  });
  return (
    <div
      className="ff-node-value-input"
      data-testid={tid('node-value-input')}
      ref={ref}
      onPointerDown={(e) => e.stopPropagation()}
      onDoubleClick={(e) => e.stopPropagation()}
    />
  );
}

/** Output nodes render only as Outputs-strip chips — the rete node view is an
 *  empty hidden placeholder that never paints or intercepts pointer events. */
function OutputRowComponent(props: NodeProps): React.JSX.Element {
  return <div style={{display: 'none'}} data-node-id={props.data.id} />;
}

interface SocketProps {
  data: TypedSocket;
}

/** Type-letter socket chip; the color travels as the `--socket-color` CSS var.
 *  Order sockets get no title so the exec-port wrapper's tooltip isn't shadowed. */
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

/** Connection painted in the source slot's color (`_color` in the payload) —
 *  the default preset hardcodes `stroke: steelblue` via styled-components. */
interface ConnectionProps {
  data: FlowScheme['Connection'] & {_color?: string; _count?: string};
}

export function FlowConnectionComponent(props: ConnectionProps): React.JSX.Element | null {
  const {path, start, end} = Presets.classic.useConnection();
  if (!path) return null;
  const color = props.data._color ?? '#8892a0';
  const count = props.data._count;

  const waypoints = (props.data as {waypoints?: Array<{x: number; y: number}>}).waypoints;
  let drawPath = path;
  if (waypoints && waypoints.length > 0 && start && end) {
    const points = [start, ...waypoints, end];
    drawPath = '';
    for (let i = 0; i < points.length - 1; i++)
      drawPath += classicConnectionPath([points[i], points[i + 1]], 0.3);
  }

  // A wide transparent "hit" path takes the pointer events; the visible 2.5px
  // path has pointerEvents:none.
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
