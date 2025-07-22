import * as bacon from 'baconjs';
import _ from 'underscore';
import {SettingsType} from './types';

/**
 * Hold on to event when holdProperty is true, and only keep them if
 * acceptProperty is true (when holdProperty becomes false).
 */
function convertToConditionalStream(valueStream: bacon.Bus<unknown>, statusStream: bacon.Bus<unknown>) {
  // Combine values with status to revert to last value when a reject is passed.
  const init: any = {
    savedValue: null,
    currentValue: null,
    lastStatus: null
  };

  const held = bacon.combineAsArray(valueStream, statusStream.toProperty(null))
    .scan(init, ({savedValue, currentValue, lastStatus}, [value, status]) => {
      // See if the status was just set
      const newStatus = lastStatus !== status;

      if (newStatus && status === 'hold') {
        // Record the currentValue as the savedValue
        return {
          savedValue: currentValue,
          currentValue,
          lastStatus: status
        };
      } else if (!newStatus && status === 'hold') {
        // Record the current value, and keep the savedValue unchanged
        return {
          savedValue,
          currentValue: value,
          lastStatus: status
        };
      } else if (newStatus && status === 'abandon') {
        // Keep the saved value
        return {
          savedValue: null,
          currentValue: savedValue,
          lastStatus: status
        };
      } else if (newStatus && status === 'accept') {
        // Keep the new value
        return {
          savedValue: null,
          currentValue,
          lastStatus: status
        };
      } else {
        // Not held, so keep the value
        return {
          savedValue: null,
          currentValue: value,
          lastStatus: status
        };
      }
    })
  // Skip the initial null value
    .skip(1)
  // Get the current value
    .map(({currentValue}) => currentValue)
  // Skip duplicate values
    .skipDuplicates()
  // property -> event stream
    .toEventStream();

  return held;
}

/**
 * Settings. A class to manage settings for a Map.
 *
 * Arguments
 * ---------
 *
 * setOption: A function, fn(key), that returns the option value for the key.
 *
 * getOption: A function, fn(key, value), that sets the option for the key and
 * value.
 *
 * conditionalOptions: The options to that are conditionally accepted when
 * changed. Changes can be abandoned by calling abandonChanges(), or accepted
 * by calling acceptChanges().
  * @param optionsWithDefaults - The current option values
  * @param conditionalOptions - The options to that are conditionally accepted
  *                             when changed. Changes can be abandoned by calling
  *                             abandon_changes(), or accepted by calling
  *                             accept_changes().
  */
export default class Settings {
  private _options: SettingsType & {[key: string]: any};
  statusBus: bacon.Bus<unknown>;
  busses: _.Dictionary<bacon.Bus<any>>;
  streams: _.Dictionary<bacon.EventStream<any>>;
  acceptedStreams: _.Dictionary<bacon.EventStream<any>>;
  get options() {
    return this._options;
  }
  constructor(optionsWithDefaults: SettingsType, conditionalOptions: (keyof SettingsType)[] = []) {
    this._options = optionsWithDefaults;

    // Manage accepting/abandoning changes
    this.statusBus = new bacon.Bus()

    // Create the options
    //@ts-ignore
    ;[this.busses, this.streams, this.acceptedStreams] = _.chain(optionsWithDefaults)
      .mapObject((value, key) => {
        const isConditional = _.contains(conditionalOptions, key);
        const {bus, stream, acceptedStream} = this.createSetting(key, value, isConditional);
        return [bus, stream, acceptedStream];
      }) // { k: [ b, s ], ... }
      .pairs() // [ [ k, [ b, s ] ], ... ]
      .map(([name, [bus, stream, acceptedStream]]) => [
        [name, bus],
        [name, stream],
        [name, acceptedStream]
      ]) // [ [ [ k, b ], [ k, s ] ], ... ]
      .unzip() // [ [ [ k, b ], ... ], [ [ k, s ], ... ] ]
      .map((x) => _.object(x)) // [ { k: b, ... }, { k: s, ... } ]
      .value();
  }

  /**
   * Set up a new bus and stream for a conditional setting (i.e. one that can be
   * canceled in the settings menu.
   */
  createSetting(name: string, initialValue: any, isConditional: boolean) {
    // Set up the bus
    const bus = new bacon.Bus();

    // Make the event stream and conditionally accept changes
    const stream = isConditional ?
      convertToConditionalStream(bus, this.statusBus) :
      bus.toEventStream();

    // Make a stream for the accepted values only. First get the latest value
    // after accepting. Also get a latest copy of the correct value if the data
    // is abandoned.
    const acceptedStream = stream.sampledBy(
      this.statusBus.filter((status) => status === 'accept' || status === 'abandon')
    ).merge(
      // Then merge with all the other changes
      stream.filter(
        this.statusBus.map((status) => status === 'accept').toProperty(true)
      )
    );

    // Get the latest
    stream.onValue((v) => { this._options[name] = v; });

    // Push the initial value
    bus.push(initialValue);

    return {bus, stream, acceptedStream};
  }

  /**
   * Deprecated. Use `set` instead.
   */
  set_conditional(name: string, value: any) { // eslint-disable-line camelcase
    console.warn('set_conditional is deprecated. Use Settings.set() instead');
    return this.set(name, value);
  }

  /**
   * Set the option. This should always be used instead of setting options
   * directly. To set options that respect the Settings menu Accept/Abandon, use
   * setConditional().
   * @param {String} name - The option name
   * @param {Any} value - The new value
   * can check whether the change was made internally to avoid loops.
   */
  set<T extends keyof SettingsType>(name: T, value: SettingsType[T]): void;
  set(name: keyof SettingsType | string, value: any): void;
  set(name: string, value: any): void {
    if (!(name in this.busses))
      throw new Error(`Invalid setting name ${name}`);

    this.busses[name].push(value);
  }

  /**
   * Deprecated. Use `get` intead.
   */
  get_option(name: string) { // eslint-disable-line camelcase
    console.warn('get_option is deprecated. Use Settings.get() instead');
    return this.get(name);
  }

  /**
   * Get an option with type inference
   * @param name - The name of the setting to get
   * @returns The value of the setting with its correct type
   */
  get<K extends keyof SettingsType>(name: K): SettingsType[K];
  get(name: string): any;
  get(name: keyof SettingsType | string): any {
    return this._options[name];
  }

  holdChanges() {
    this.statusBus.push('hold');
  }

  abandonChanges() {
    this.statusBus.push('abandon');
  }

  acceptChanges() {
    this.statusBus.push('accept');
  }
}
