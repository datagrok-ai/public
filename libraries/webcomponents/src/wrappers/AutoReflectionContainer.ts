import {Viewer} from 'datagrok-api/dg';
import {InstancesContainer} from './InstancesContainer';
import {map} from 'rxjs/operators';

// Uses Viewer dart reflection to automatically wire viewer instances
// to a webcomponet container.
export abstract class AutoReflectionContainer<T = any> extends InstancesContainer<Viewer<T>> {
  constructor() {
    super();
  }

  override applyInstanceProps(instance: Viewer<T>, props: Partial<T>) {
    instance.setOptions(props);
  }

  override getInstancePropertiesDesriptions(instance: Viewer<T>) {
    const dartProps = instance.getProperties();
    const props = Object.fromEntries(dartProps.map((prop) => {
      const name = prop.name;
      const value = prop.propertyType;
      return [name, value];
    }));
    return {
      type: instance.type,
      props,
    };
  }

  override getInstanceProperties(instance: Viewer<T>) {
    const jsOptions = instance.getOptions(true)?.look;
    const dartProps = instance.getProperties();
    return Object.fromEntries(dartProps.map((prop) => {
      const name = prop.name;
      const value = jsOptions[name];
      return [name, value];
    }));
  }

  override getInstanceProperty(instance: Viewer<T>, name: any) {
    const jsOptions = instance.getOptions(true)?.look;
    return jsOptions[name];
  }

  override getInstanceEvents(instance: Viewer<T>) {
    return instance.onDartPropertyChanged.pipe(
      map((dev: any) => {
        const name = dev.dart.c;
        return {
          name,
          data: this.getInstanceProperty(instance, name),
        };
      }),
    );
  }
}
