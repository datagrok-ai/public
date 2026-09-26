/* eslint-disable camelcase */
/** @jsx h */
import {h, Component} from 'preact';
import {Range} from 'rc-slider-preact';
import 'rc-slider-preact/lib/index.css';
import './css/DistributionTooltipComponent.css';
import * as utils from './ts/utils';
import { CobraModelData, FluxHistogram, ReactionSamplingDistribution } from './ts/types';

// const utils = require('escher/src/utils.js')
const WIDTH = 320;
const HEIGHT = 175;
const NON_DIST_HEIGHT = 175;
// or: import { WIDTH } from './constants'

const tooltipStyle = {
  display: '-webkit-flex',
  flexDirection: 'column',
  boxSizing: 'border-box',
  width: WIDTH + 'px',
  borderRadius: '2px',
  border: '1px solid #b58787',
  padding: '1.5%',
  backgroundColor: 'rgba(255, 255, 255, 1)',
  textAlign: 'left',
  fontSize: '16px',
  fontFamily: 'sans-serif',
  color: '#111',
  boxShadow: '4px 6px 20px 0px rgba(0, 0, 0, 0.4)',
  position: 'relative',
  zIndex: '3'
};
const fluxDisplayStyle = {
  alignSelf: 'center',
  color: 'black',
  fontSize: '12px',
  fontWeight: 'bold',
  visibility: 'visible'
};
const indicatorStyle = {
  marginTop: '0.5%',
  display: '-webkit-flex',
  flexDirection: 'column',
  alignSelf: 'center',
  width: '14px',
  height: '14px',
  visibility: 'visible'
};

/** Axis label with just enough decimals to tell values `span` apart (6 for a single value). */
function formatFlux(value: number, span: number): string {
  const decimals = span > 0 ? Math.min(8, Math.max(0, Math.ceil(-Math.log10(span)) + 2)) : 6;
  return String(parseFloat(value.toFixed(decimals)));
}

/** Histogram of one reaction's sampled fluxes, over that reaction's own [min, max]. */
export class HistogramComponent extends Component<{histogram: FluxHistogram}> {
  render() {
    const {min, max, counts} = this.props.histogram;
    const height = 150; // Height in pixels
    const padding = 20; // Padding for the axis labels
    const span = max - min;
    // every sample had the same flux: one centered bar, labeled with that value
    const flat = !(span > 0);
    const barWidth = flat ? 6 : 100 / counts.length;
    const binWidth = flat ? 0 : span / counts.length;
    const ticks = flat ? [min] : [min, (min + max) / 2, max];

    const maxCount = counts.reduce((m, c) => Math.max(m, c), 0);
    const scale = maxCount > 0 ? (height - padding * 2) / maxCount : 0;

    return (
      <div className='escher-distribution-histogram' style={{
        position: 'relative',
        width: '100%',
        height: `${height}px`,
        paddingTop: `${padding}px`,

      }}>
        {/* Render bars */}
        <div style={{display: 'flex', justifyContent: 'center', height: `${height - padding * 2}px`, alignItems: 'flex-end', width: '100%'}}>
          {counts.map((count, index) => (
            <div
              key={index}
              title={(flat ? formatFlux(min, 0) :
                `${formatFlux(min + index * binWidth, binWidth)} to ${formatFlux(min + (index + 1) * binWidth, binWidth)}`) +
                `: ${count} sample${count === 1 ? '' : 's'}`}
              style={{
                width: `${barWidth}%`,
                height: `${count * scale}px`,
                backgroundColor: '#2083D590',
                marginRight: index === counts.length - 1 ? '0' : '1px',
                borderBottom: '1px solid black'
              }}
            />
          ))}
        </div>

        {/* X-axis labels */}
        <div style={{
          display: 'flex',
          justifyContent: flat ? 'center' : 'space-between',
          width: '100%',
          position: 'absolute',
          left: '0'
        }}>
          {ticks.map((t) => <div>{formatFlux(t, span)}</div>)}
        </div>
      </div>
    );
  }
}

export type DistributionTooltipProps = {
  type: string;
  name: string;
  model: CobraModelData;
  oldModel: CobraModelData;
  biggId: string;
  reactionData: { [key: string]: number };
  upperRange: number;
  lowerRange: number;
  step: number;
  fluxDistributions: ReactionSamplingDistribution;
  samplingFunction: () => void;
  objectives: { [key: string]: number };
  sliderChange: (bounds: number[], biggId: string) => void;
}

class DistributionTooltipComponent extends Component<DistributionTooltipProps> {
  state: {
    lowerBound: number;
    upperBound: number;
    lowerBoundOld: number;
    upperBoundOld: number;
    currentFlux: number;
    coefficient: number;
    lowerBoundString: string;
    upperBoundString: string;
    type: string;
    name: string;
    reactionInModel: boolean;
    indicatorStyle: { [key: string]: string };
    fluxDisplayStyle: { [key: string]: string };
    tooltipStyle: { [key: string]: string };
    objectives: { [key: string]: number };
    fluxDistributions: ReactionSamplingDistribution;
  };

  constructor(props: DistributionTooltipProps) {
    super(props);
    this.state = {
      lowerBound: 0,
      upperBound: 0,
      lowerBoundOld: 0,
      upperBoundOld: 0,
      currentFlux: 0,
      coefficient: 0,
      lowerBoundString: '',
      upperBoundString: '',
      type: 'reaction',
      name: '',
      reactionInModel: true,
      indicatorStyle,
      fluxDisplayStyle,
      tooltipStyle,
      objectives: {},
      fluxDistributions: props.fluxDistributions // should be a Map
    };
    this.initState(props, {});
  }

  componentWillReceiveProps(nextProps) {
    this.initState(nextProps, this.props);
  }

  initState(nextProps: DistributionTooltipProps, lastProps: Partial<DistributionTooltipProps>) {
    // By default, reaction is not in model
    let reactionInModel = false;

    // If the selected map object is a reaction and there is a model, collect the necessary
    // flux data and calculate placement of arrows and current flux label.
    // Otherwise, only pass the type to the tooltip.
    if (nextProps.type === 'reaction' &&
        !(nextProps.model === undefined || nextProps.model == null)) {
      const fluxData: Partial<typeof this.state> = {};

      // Only updates all of the flux data when the reaction, model, or objective changes.
      // Otherwise only updates the current flux.
      if (nextProps.biggId !== lastProps.biggId ||
        nextProps.model !== lastProps.model) {
        //
        for (let i = 0, l = nextProps.model.reactions.length; i < l; i++) {
          //
          if (nextProps.model.reactions[i].id === nextProps.biggId) {
            reactionInModel = true;
            fluxData.lowerBound = nextProps.model.reactions[i].lower_bound;
            fluxData.upperBound = nextProps.model.reactions[i].upper_bound;
            fluxData.lowerBoundString = nextProps.model.reactions[i].lower_bound.toString();
            fluxData.upperBoundString = nextProps.model.reactions[i].upper_bound.toString();
            fluxData.lowerBoundOld = nextProps.oldModel.reactions[i].lower_bound;
            fluxData.upperBoundOld = nextProps.oldModel.reactions[i].upper_bound;
            fluxData.name = nextProps.model.reactions[i].name;
            fluxData.coefficient = nextProps.model.reactions[i].objective_coefficient;
            if (nextProps.reactionData != null)
              fluxData.currentFlux = nextProps.reactionData[nextProps.biggId];
            else
              fluxData.currentFlux = null;

            break;
          }
        }
      } else {
        reactionInModel = true;
        for (let i = 0, l = nextProps.model.reactions.length; i < l; i++) {
          if (nextProps.model.reactions[i].id === nextProps.biggId) {
            if (nextProps.reactionData != null)
              fluxData.currentFlux = nextProps.reactionData[nextProps.biggId];

            break;
          }
        }
      }
      // For calculating placement of the current flux indicator arrow and label
      let textOffset = {};
      let arrowPosition = {};
      if (nextProps.reactionData != null) {
        if (Math.abs(fluxData.currentFlux) > nextProps.upperRange) {
          arrowPosition = {
            marginLeft: fluxData.currentFlux / Math.abs(fluxData.currentFlux) * 100 + '%'
          };
        } else {
          arrowPosition = {
            marginLeft: fluxData.currentFlux / (nextProps.upperRange + 1) * 100 + '%'
          };
        }
        if (fluxData.currentFlux > 0.65 * nextProps.upperRange)
          textOffset = {alignSelf: 'flex-end'};
        else if (fluxData.currentFlux < 0.65 * nextProps.lowerRange)
          textOffset = {alignSelf: 'flex-start'};
        else
          textOffset = {marginLeft: fluxData.currentFlux / (nextProps.upperRange + 1) * 100 + '%'};
      } else {
        textOffset = {
          visibility: 'hidden'
        };
        arrowPosition = {
          visibility: 'hidden'
        };
      }
      this.setState({
        ...fluxData,
        reactionInModel,
        fluxDisplayStyle: {...fluxDisplayStyle, ...textOffset},
        indicatorStyle: {...indicatorStyle, ...arrowPosition},
        type: nextProps.type,
        objectives: {}
      });
    } else {
      this.setState({
        type: nextProps.type,
        reactionInModel: false
      });
    }
  }

  get_size() {
    return {width: WIDTH, height: HEIGHT};
  }

  /**
   * Due to bugs inherent in the slider, all values must be converted to fall
   * somewhere onto the positive number line in order to be displayed correctly
   * on the slider. This function takes a given value and converts it so that it
   * will display on the number line correctly.
   * @param {number} value - Physiologically relevant number to be displayed on
   * slider.
   */
  fluxConverter(value: number) {
    return value < this.props.lowerRange ? // Add parenthesis for better readability
      -1000 :
      value > this.props.upperRange ?
        1000 :
        value + this.props.upperRange + 1;
  }

  /**
   * Due to bugs inherent in the slider, all values must be converted to fall
   * somewhere onto the positive number line in order to be displayed correctly
   * on the slider. This function takes values from the slider and converts them
   * back to a physiologically relevant value.
   * @param {number} value - Slider value to be converted to physiologically
   * relevant value.
   */
  tipConverter(value: number): number {
    let sigFig = 0;
    if (this.props.step < 1)
      sigFig = Math.ceil(-Math.log10(this.props.step));

    return value < 1 ?
      -1000 :
      value > (2 * this.props.upperRange + 1) ?
        1000 :
        parseFloat((value - (this.props.upperRange + 1)).toFixed(sigFig));
  }

  /**
   * Function for applying tipConverter to arrays of numbers.
   * @param {number[]} array - Pair of values (lower and upper bounds,
   * respectively) to be converted from slider values to physiologically
   * relevant values.
   */
  boundConverter(array: (number | string)[]): number[] {
    return array.map(this.tipConverter.bind(this));
  }

  onTouchStart() {
    const dragStyle = {...this.state.tooltipStyle, backgroundColor: 'rgba(255, 255, 255, 0.8)'};
    this.setState({
      tooltipStyle: dragStyle
    });
  }

  handleMarkerPosition(currentFlux: number) {
    if (this.props.reactionData != null) {
      return currentFlux < this.props.lowerRange ? // Add parenthesis for better readability
        0 :
        currentFlux > this.props.upperRange ?
          2 * (this.props.upperRange + 1) :
          currentFlux + this.props.upperRange + 1;
    } else
      return this.props.upperRange + 1;
  }

  handleKeyUp(event: Event, bounds: string[]) {
    this.setState({
      lowerBoundString: bounds[0],
      upperBoundString: bounds[1]
    });
    if (parseFloat(bounds[0]) !== this.state.lowerBound || parseFloat(bounds[1]) !== this.state.upperBound) {
      if (isNaN(parseFloat((event.target as HTMLInputElement).value)))
        console.error('Invalid Bounds');
      else
        this.sliderChange([parseFloat(this.state.lowerBoundString), parseFloat(this.state.upperBoundString)]);
    }
  }

  sliderChange(bounds: number[]) {
    if (bounds[0] !== parseFloat(this.state.lowerBoundString) || bounds[1] !== parseFloat(this.state.upperBoundString)) {
      this.setState({
        lowerBound: bounds[0],
        upperBound: bounds[1],
        lowerBoundString: bounds[0].toString(),
        upperBoundString: bounds[1].toString()
      });
    } else {
      this.setState({
        lowerBound: bounds[0],
        upperBound: bounds[1]
      });
    }
    this.props.sliderChange(bounds, this.props.biggId);
  }

  onTouchEnd(bounds: number[]) {
    this.sliderChange(bounds);
    this.setState({
      tooltipStyle
    });
  }

  resetReaction() {
    this.sliderChange([this.state.lowerBoundOld, this.state.upperBoundOld]);
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
    const url = (type === 'gene' ?
      pref + 'search?query=' + biggId :
      pref + 'universal/' + type + 's/' + this.decompartmentalizeCheck(biggId, type));
    window.open(url);
  }

  render() {
    if (this.state.type === 'reaction' && this.state.reactionInModel) {
      // get state of max and min buttons
      const minimizeActive = this.props.objectives[this.props.biggId] === -1;
      const minimizeDisabled = Object.keys(this.props.objectives).length === 1 && minimizeActive;
      const maximizeActive = this.props.objectives[this.props.biggId] === 1;
      const maximizeDisabled = Object.keys(this.props.objectives).length === 1 && maximizeActive;
      const fluxHistogram = this.props.fluxDistributions?.data?.get(this.props.biggId);

      // create a reaction string from the model
      const reaction = this.props.model?.reactions.find((r) => r.id === this.props.biggId);
      let reactionString: string | null = null;
      if (reaction && reaction.metabolites) {
        const inputs = (Object.entries(reaction.metabolites).filter(([_, coeff]) => coeff < 0)).map(([metId, coeff]) => Math.abs(coeff) !== 1 ? `${Math.abs(coeff)}*${metId}` : metId);
        const outputs = (Object.entries(reaction.metabolites).filter(([_, coeff]) => coeff > 0)).map(([metId, coeff]) => Math.abs(coeff) !== 1 ? `${Math.abs(coeff)}*${metId}` : metId);
        const isReversible = reaction.reversibility || (reaction.upper_bound > 0 && reaction.lower_bound < 0);

        reactionString = `${inputs.join(' + ')} ${isReversible ? '⇋' : '→'} ${outputs.join(' + ')}`;
        
      }

      return (
        <div className='escher-tooltip'
          style={{
            ...this.state.tooltipStyle,
            touchAction: 'none'
          }}
        >
          <div className='biggId'>
            {this.props.biggId}
          </div>
          <div className='name'>
            {this.props.name}
          </div>
          {reactionString && (
            <div className='reactionString'>
              {reactionString}
            </div>
          )}
          <div className='slider' style={{width: WIDTH - 22 + 'px'}}>
            <Range
              onBeforeChange={() => this.onTouchStart()}
              style={{alignSelf: 'center'}}
              min={0}
              max={2 * (this.props.upperRange + 1)}
              step={this.props.step}
              value={[
                this.fluxConverter(this.state.lowerBound),
                this.fluxConverter(this.state.upperBound)
              ]}
              // marks={{[this.fluxConverter(this.state.currentFlux)]: ''}} // no current flux
              tipFormatter={() => null}
              allowCross={false}
              pushable={0}
              onChange={(f) => this.sliderChange(this.boundConverter(f as number[]))}
              onAfterChange={(f) => this.onTouchEnd(this.boundConverter(f as number[]))}
            />
            {this.state.currentFlux != null && (
              <div style={{width: '100%', position: 'relative', display: 'flex', flexDirection: 'column'}}>
                <div className='indicator' style={this.state.indicatorStyle}>
                <svg viewBox='0 0 100 100' height='100%' width='100%'>
                  <defs>
                    <marker id='markerArrow1' viewBox='0 0 6 6' refX='4' refY='3' orient='auto'>
                      <path d='M5,3 L3,5 L3,1 Z' fill='black' stroke='black' />
                    </marker>
                  </defs>
                  <line x1='50' y1='75' x2='50' y2='20' stroke-width='25' stroke='black' marker-end={'url(#markerArrow1)'} />
                </svg>
              </div>
              <div className='fluxDisplay' style={this.state.fluxDisplayStyle}>
                Flux: {this.state.currentFlux!.toFixed(2)}
              </div>
            </div>
            )}
            
          </div>
          {/* Kebab case for class names?  */}

          <div className='interfacePanel'>
            <div className='labels'>
              <div
                style={{
                  //  float: 'left',
                  fontSize: '12px'
                }}
              >
              Lower bound
              </div>
              <div
                style={{
                  //  float: 'right',
                  fontSize: '12px'
                }}
              >
              Upper bound
              </div>
            </div>


            <div className='inputPanel'>
              <input
                type='text'
                className='input'
                value={this.state.lowerBoundString}
                onFocus={(event) => (event.target as HTMLInputElement).select()}
                onKeyUp={
                  (event) => this.handleKeyUp(event, [(event.target as HTMLInputElement).value, this.state.upperBoundString])
                }
              />
              <input
                type='text'
                className='input'
                value={this.state.upperBoundString}
                onFocus={(event) => (event.target as HTMLInputElement).select()}
                onKeyUp={
                  (event) => this.handleKeyUp(event, [this.state.lowerBound.toString(), (event.target as HTMLInputElement).value?.toString()])
                }
              />
            </div>
            <div className='buttonBar'>
              <button
                className='button'
                onClick={
                  () => this.sliderChange([0, 0])
                }
              >
                Knockout
              </button>
              <button
                className='button'
                onClick={() => this.resetReaction()}
              >
                Reset
              </button>

              {this.props.samplingFunction != null && (
                <button
                  className='button'
                  onClick={() => this.props.samplingFunction()}
                >
                  Run Sampling
                </button>
              )

              }

              {/* <button
                className={maximizeActive ? 'active' : ''}
                onClick={() => this.props.setObjective(this.props.biggId, 1)}
                disabled={maximizeDisabled}
              >
                Maximize
              </button>
              <button
                className={minimizeActive ? 'active' : ''}
                onClick={() => this.props.setObjective(this.props.biggId, -1)}
                disabled={minimizeDisabled}
              >
                Minimize
              </button> */}
            </div>
          </div>
          { fluxHistogram && (
            <div className="histogramPanel">
              <HistogramComponent histogram={fluxHistogram} />
            </div>
          )}
        </div>
      );
    } else if (this.state.type === 'reaction' && !this.state.reactionInModel) {
      return (
        <div className='escher-tooltip'
          style={{
            ...tooltipStyle,
            height: '50px',
            width: '400px'
          }}
        >
          <div className='biggId'>
            {this.props.biggId} is not in the model
          </div>
        </div>
      );
    } else {
      const decomp = this.decompartmentalizeCheck(this.props.biggId, this.props.type);
      return (
        <div
          className='escher-tooltip'
          style={{
            ...this.state.tooltipStyle,
            touchAction: 'none',
            height: 'fit-content',
            width: '250px'
          }}
        >
          <div className='biggId'>
            {this.props.biggId}
          </div>
          <div className='name'>
            {this.props.name}
          </div>
          <button
            className='button'
            style={{width: 'fit-content', marginTop: '10px'}}
            onClick={() => this.openBigg()}
          >
            Open {decomp} in BiGG Models
          </button>
        </div>
      );
    }
  }
}

export default DistributionTooltipComponent;
