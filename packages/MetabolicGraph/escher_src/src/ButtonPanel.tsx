/** @jsx h */
import {h, Component} from 'preact';
import './css/ButtonPanel.css';
import Settings from './ts/Settings';
import {EscherMap} from './ts/escherMap';
import ZoomContainer from './ts/ZoomContainer';
import BuildInput from './ts/BuildInput';

/**
 * ButtonPanel. Sets up the button panel for Builder. Currently calls a
 * re-render on itself every time the mode is changed. This can be removed upon
 * porting Builder to Preact
 */
class ButtonPanel extends Component<{
  settings: Settings;
  map: EscherMap;
  zoomContainer: ZoomContainer;
  full_screen: () => void;
  mode: string;
  isFullScreen: boolean;
  setMode: (mode: string) => void;
  buildInput: BuildInput;
}> {
  render() {
    const menuSetting = this.props.settings.get('menu');
    const enableKeys = this.props.settings.get('enable_keys');
    const enableEditing = this.props.settings.get('enable_editing');

    return (
      <ul className='button-panel'>
        {!this.props.settings.get('pathFindingDisabled') && (
          <li>
            <button
              className='button btn'
              onClick={() => this.props.map.findAndHighlightShortest()}
              title={'Find Shortest path'}
            >
              <i className='icon-shortest-path' />
            </button>
            <div className='kth-reaction-arrows'>
              <button
                className='buttonGroup btn'
                title={`Previous shortest reaction`}
                onClick={() => this.props.map.findAndHighlightShortest(-1)}
              >
                <i className='icon-left-big' />
              </button>

              <button
                className='buttonGroup btn'
                title={`Next shortest reaction`}
                onClick={() => this.props.map.findAndHighlightShortest(1)}
              >
                <i className='icon-right-big' />
              </button>

            </div>

          </li>)}

        {!this.props.settings.get('pathFindingDisabled') && (
          <li>
            <button
              className='button btn'
              onClick={() => this.props.map.findKthOutOrIngoingReactions(0, true)}
              title={'Highlight outgoing reactions'}
            >
              <i className='icon-arrows-outside' />
            </button>
            <div className='kth-reaction-arrows'>
              <button
                className='buttonGroup btn'
                title={`Decrement outgoing reactions`}
                onClick={() => this.props.map.findKthOutOrIngoingReactions(-1, true)}
              >
                <i className='icon-left-big' />
              </button>

              <button
                className='buttonGroup btn'
                title={`Increment outgoing reactions`}
                onClick={() => this.props.map.findKthOutOrIngoingReactions(1, true)}
              >
                <i className='icon-right-big' />
              </button>

            </div>

          </li>)}

        {!this.props.settings.get('pathFindingDisabled') && (
          <li>
            <button
              className='button btn'
              onClick={() => this.props.map.findKthOutOrIngoingReactions(0, false)}
              title={'Highlight ingoing reactions'}
            >
              <i className='icon-arrows-inside' />
            </button>
            <div className='kth-reaction-arrows'>
              <button
                className='buttonGroup btn'
                title={`Decrement ingoing reactions`}
                onClick={() => this.props.map.findKthOutOrIngoingReactions(-1, false)}
              >
                <i className='icon-left-big' />
              </button>

              <button
                className='buttonGroup btn'
                title={`Increment ingoing reactions`}
                onClick={() => this.props.map.findKthOutOrIngoingReactions(1, false)}
              >
                <i className='icon-right-big' />
              </button>

            </div>

          </li>)}

        <li>
          <button
            className='button btn'
            onClick={() => this.props.zoomContainer.zoom_in()}
            title={`Zoom in${enableKeys ? ' (+)' : ''}`}
          >
            <i className='icon-zoom-in' />
          </button>
        </li>
        <li>
          <button
            className='button btn'
            onClick={() => this.props.zoomContainer.zoom_out()}
            title={`Zoom out${enableKeys ? ' (-)' : ''}`}
          >
            <i className='icon-zoom-out' />
          </button>
        </li>
        <li>
          <button
            className='button btn'
            onClick={() => this.props.map.zoom_extent_canvas()}
            title={`Zoom to canvas${enableKeys ? ' (1)' : ''}`}
          >
            <i className='icon-resize-full' />
          </button>
        </li>
        <li style={{display: this.props.settings.get('full_screen_button') !== false ? 'block' : 'none'}}>
          <button
            className={`button btn ${this.props.isFullScreen ? 'active-button' : ''}`}
            onClick={() => this.props.full_screen()}
            title={'Toggle full screen'}
          >
            <i className='icon-resize-full-alt' />
          </button>
        </li>
        <li
          className='grouping'
          style={{display: menuSetting === 'all' && enableEditing ? 'block' : 'none'}}
        >
          <button
            className='buttonGroup btn'
            title={`Pan mode${enableKeys ? ' (Z)' : ''}`}
            for='zoom'
            id={this.props.mode === 'zoom' ? 'currentMode' : null}
            onClick={() => this.props.setMode('zoom')}
          >
            <i className='icon-move' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Select mode${enableKeys ? ' (V)' : ''}`}
            for='brush'
            id={this.props.mode === 'brush' ? 'currentMode' : null}
            onClick={() => this.props.setMode('brush')}
          >
            <i className='icon-mouse-pointer' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Add reaction mode${enableKeys ? ' (N)' : ''}`}
            for='build'
            onClick={() => this.props.setMode('build')}
            id={this.props.mode === 'build' ? 'currentMode' : null}>
            <i className='icon-wrench' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Rotate mode${enableKeys ? ' (R)' : ''}`}
            for='rotate'
            id={this.props.mode === 'rotate' ? 'currentMode' : null}
            onClick={() => this.props.setMode('rotate')}
          >
            <i className='icon-cw' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Text mode${enableKeys ? ' (T)' : ''}`}
            for='text'
            id={this.props.mode === 'text' ? 'currentMode' : null}
            onClick={() => this.props.setMode('text')}
          >
            <i className='icon-font' />
          </button>
        </li>
        <li
          className='grouping'
          style={{display: this.props.mode === 'build' && menuSetting === 'all' && enableEditing ? 'block' : 'none'}}
        >
          <button
            className='buttonGroup btn'
            title={`Direction arrow${enableKeys ? ' (←)' : ''}`}
            onClick={() => this.props.buildInput.direction_arrow.left()}
          >
            <i className='icon-left-big' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Direction arrow${enableKeys ? ' (→)' : ''}`}
            onClick={() => this.props.buildInput.direction_arrow.right()}
          >
            <i className='icon-right-big' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Direction arrow${enableKeys ? ' (↑)' : ''}`}
            onClick={() => this.props.buildInput.direction_arrow.up()}
          >
            <i className='icon-up-big' />
          </button>
          <button
            className='buttonGroup btn'
            title={`Direction arrow${enableKeys ? ' (↓)' : ''}`}
            onClick={() => this.props.buildInput.direction_arrow.down()}
          >
            <i className='icon-down-big' />
          </button>
        </li>
      </ul>
    );
  }
}

export default ButtonPanel;
