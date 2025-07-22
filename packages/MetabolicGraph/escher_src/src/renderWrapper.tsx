/** @jsx h */
import {h, Component, render, Ref} from 'preact';

class Wrapper extends Component {
  state: {display?: boolean} = {display: false};
  props: {component: any, connectSetStateFn: (props: any) => void, ref: any, refPassthrough: any};
  constructor(props: {component: any, connectSetStateFn: (props: any) => void, ref: any, refPassthrough: any}) {
    super(props);
  }

  componentDidMount() {
    this.props.connectSetStateFn((props) => this.setState(props));
  }

  is_visible() { // eslint-disable-line camelcase
    return this.state.display;
  }

  render() {
    if (!this.state.display) return null;

    // Pass the new props, and always pass the ref
    return (
      <this.props.component
        setDisplay={(display) => this.setState({display})}
        ref={this.props.refPassthrough}
        {...this.state}
      />
    );
  }
}

/**
 * Wrapper for better integration of Preact components with Escher. The
 * component can be updated using the connectSetStateFn to set up a callback for
 * updates from other components.
 * @param {} component - A Preact component
 *
 * @param {} ref - A preact ref for the wrapper so that the "display" state can
 *                 be tracked.
 */
function renderWrapper(
  component: any,
  ref: any,
  connectSetStateFn: (props: any) => void,
  divNode: Element,
  refPassthrough: Ref<any> | null = null
) {
  render(
    <Wrapper
      component={component}
      connectSetStateFn={connectSetStateFn}
      ref={ref}
      refPassthrough={refPassthrough}
    />,
    divNode,
    // If there is already a div, re-render it. Otherwise make a new one
    divNode.children.length > 0 ? divNode.firstChild as Element : undefined
  );
}

export default renderWrapper;
