/* global window */
/** @jsx h */

/**
 * Define a Tooltip component and interface with Preact.
 */
import {h, Component} from 'preact';
import './css/DefaultTooltip.css';
import * as utils from './ts/utils';

class DefaultTooltip extends Component<{
  biggId: string;
  name: string;
  data: string;
  type: string;
  disableTooltips: () => void;
}> {
  constructor() {
    super();
    this.openBigg = this.openBigg.bind(this);
  }

  decompartmentalizeCheck(id: string, type: string) {
  // ID without compartment, if metabolite.
    return type === 'metabolite' ?
      utils.decompartmentalize(id)[0] :
      id;
  }

  openBigg() {
    const type = this.props.type;
    const biggId = this.props.biggId;
    const pref = 'http://bigg.ucsd.edu/';
    const url = type === 'gene' ?
      `${pref}search?query=${biggId}` :
      `${pref}universal/${type}s/${this.decompartmentalizeCheck(biggId, type)}`;
    window.open(url);
  }

  capitalizeFirstLetter(s: string) {
    if (typeof s !== 'string') {
      console.warn('capitalizeFirstLetter was passed something other than a string');
      return '';
    }
    return (s.charAt(0).toUpperCase() + s.slice(1));
  }

  render() {
    const decomp = this.decompartmentalizeCheck(this.props.biggId, this.props.type);
    const biggButtonText = `Open ${decomp} in BiGG Models.`;
    return (
      <div className='default-tooltip'>
        <div className='id'>
          {this.props.biggId}
        </div>
        <div className='name'>
          name: {this.props.name}
        </div>
        <div className='data'>
          data: {
            this.props.data && this.props.data !== '(nd)' ? this.props.data : 'no data'
          }
        </div>
        <button onClick={this.openBigg}>
          {biggButtonText}
        </button>
        <div className='top-right'>
          <div className='type-label'>
            {this.capitalizeFirstLetter(this.props.type)}
          </div>
          <a onClick={this.props.disableTooltips}>
            Disable Tooltips
          </a>
        </div>
      </div>
    );
  }
}

export default DefaultTooltip;
