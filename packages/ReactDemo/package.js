class ReactDemoPackage extends DG.Package {

    //name: ReactInfoPanel
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    testInfoPanel(smiles) {
        return ui.Widget.react(React.createElement(Hello, {toWhat: smiles}, null));
    }

    //name: TestDialog
    testDialog() {
        var reactHost = ui.div();
        ReactDOM.render(React.createElement(Hello, {toWhat: 'World'}, null), reactHost);

        ui.dialog('React custom components')
            .add(reactHost)
            .show();
    }
}


class Hello extends React.Component {
    render() {
        return React.createElement('div', null, `Hello ${this.props.toWhat}`);
    }
}