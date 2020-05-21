// Different ways to render a molecule using the openchemlib library
//
// For more details, visit https://github.com/cheminfo/openchemlib-js

grok.newView('chem rendering').root.appendChild(ui.div([
    ui.h1('SVG rendering'),
    grok.chem.svgMol('c1(ccc2N=C(C)N(C(=O)c2c1)c3ccc(OC)cc3)NC(=S)Nc4ccccc4'),
    grok.chem.svgMol('c1(ccc2N=C(C)N(C(=O)c2c1)c3ccc(OC)cc3)NC(=S)Nc4ccccc4', 400, 300)
]));
