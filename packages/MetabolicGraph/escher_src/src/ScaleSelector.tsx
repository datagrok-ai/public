/** @jsx h */

import {h, Component} from 'preact';

class ScaleSelector extends Component<{
  disabled: boolean;
  children: any;
}> {
  state: {
    visible: boolean;
  };
  wrapperRef: HTMLElement | null = null;

  constructor(props: {disabled: boolean; children: any}) {
    super(props);
    this.setWrapperRef = this.setWrapperRef.bind(this);
    this.handleClickOutside = this.handleClickOutside.bind(this);
    this.state = {
      visible: false
    };
  }

  componentDidMount() {
    document.addEventListener('mouseup', this.handleClickOutside);
  }

  // Reference for hiding the menu when a mouse event happens outside
  setWrapperRef(node: HTMLElement) {
    this.wrapperRef = node;
  }

  handleClickOutside(event: MouseEvent) {
    if (this.wrapperRef && !this.wrapperRef.contains(event.target as Node))
      this.setState({visible: false});
  }

  render() {
    return (
      <div className='selector'>
        <div
          className={
            [
              'selectorTitle',
              (this.props.disabled ? 'disabled' : '')
            ].join(' ')
          }
          ref={this.setWrapperRef}
          onClick={() => {
            if (!this.props.disabled)
              this.setState({visible: !this.state.visible});
          }}
        >
          Preset Scale Selections
          <i className='icon-sort-down' />
        </div>
        <div
          className='selectorMenu'
          style={
            this.state.visible ?
              {display: 'block'} :
              {display: 'none'}}
        >
          {this.props.children.map((listItem: any) => {
            return listItem;
          })}
        </div>
      </div>
    );
  }
}

export default ScaleSelector;
