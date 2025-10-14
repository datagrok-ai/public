/**
 * MenuButton. Handles the individual menu items within the Dropdown menus.
 * Takes in a name and a function and binds the function to a button. If the
 * type prop is defined as 'load' will instead render a label with an attached
 * hidden input[file] element.
 */

/** @jsx h */
import {h, Component} from 'preact';
import _ from 'underscore';
import * as utils from './ts/utils';
import * as dataStyles from './ts/dataStyles';

class MenuButton extends Component<
{
  name: string;
  onClick: (data: any) => void;
  assignKey?: (fn: () => void) => void;
  type?: string;
  disabledButtons?: string[];
  checkMark?: boolean;
}> {
  constructor(props) {
    super(props);
    this.assignKeyForInput = this.assignKeyForInput.bind(this);
  }

  handleFileInput(target: HTMLInputElement) {
    const file = target.files[0];
    const reader = new window.FileReader();
    reader.onload = () => {
      utils.load_json_or_csv(file, dataStyles.csv_converter, (e, d) => this.props.onClick(d));
    };
    if (file !== undefined)
      reader.readAsText(file);

    // reset input
    target.value = null;
  }

  assignKeyForInput(node) {
    if (this.props.assignKey)
      this.props.assignKey(() => node.click());
  }

  render() {
    const disabled = _.contains(this.props.disabledButtons, this.props.name.replace(/ \(.*\)$/, ''));
    if (this.props.type === 'load') {
      return (
        <label
          className='menuButton'
          tabIndex={disabled ? -1 : 0}
          id={disabled ? 'disabled' : ''}
        >
          <input
            type='file'
            onChange={(event) => this.handleFileInput(event.target as HTMLInputElement)}
            disabled={disabled}
            ref={this.assignKeyForInput}
          />
          {this.props.name}
        </label>
      );
    } else if (this.props.checkMark) {
      return (
        <li
          className='menuButton'
          tabIndex={disabled ? -1 : 0}
          onClick={this.props.onClick}
          id={disabled ? 'disabled' : ''}
        >
          <i className='icon-ok' aria-hidden='true'>&nbsp;</i>
          {this.props.name}
        </li>
      );
    } else {
      return (
        <li
          className='menuButton'
          tabIndex={disabled ? -1 : 0}
          onClick={disabled ? null : this.props.onClick}
          id={disabled ? 'disabled' : ''}
        >
          {this.props.name}
        </li>
      );
    }
  }
}

export default MenuButton;
