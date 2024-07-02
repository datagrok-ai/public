/*******************************************************************************
* Copyright (C)2018, The Pistoia Alliance
*  Version 1.1.0.2018--03-13
* 
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

// https://github.com/PistoiaHELM/HELMEditor/blob/master/resources/conf/DefaultMonomerCategorizationTemplate.xml
// 

/**
@project HELM Web Editor
@version 1.1.2
@description HELM Web Editor built on JSDraw.Lite
*/

/**
* HELM namespace
* @namespace org.helm.webeditor
*/

/**
* HELM Version
* @property org.helm.webeditor.version
*/


if (typeof (org) == "undefined")
    org = {};
if (org.helm == null)
    org.helm = {};

org.helm.webeditor = {
    kVersion: "1.1.2",
    atomscale: 2,
    bondscale: 1.6,
    allowHeadToHeadConnection: true,
    ambiguity: false,

    HELM: {
        BASE: "HELM_BASE",
        SUGAR: "HELM_SUGAR",
        LINKER: "HELM_LINKER",
        AA: "HELM_AA",
        CHEM: "HELM_CHEM",
        BLOB: "HELM_BLOB",
        NUCLEOTIDE: "HELM_NUCLETIDE" // only for the combo *
    },

    blobtypes: ["Bead", "Gold Particle"],

    /**
    * Test if a node is HELM monomer
    * @function isHelmNode
    */
    isHelmNode: function (a) {
        if (a == null)
            return false;

        var biotype = typeof (a) == "string" ? a : a.biotype();
        return biotype == org.helm.webeditor.HELM.BASE || biotype == org.helm.webeditor.HELM.SUGAR || biotype == org.helm.webeditor.HELM.LINKER ||
            biotype == org.helm.webeditor.HELM.AA || biotype == org.helm.webeditor.HELM.CHEM || biotype == org.helm.webeditor.HELM.BLOB || biotype == org.helm.webeditor.HELM.NUCLEOTIDE;
    },

    /**
    * List HELM Monomer Types
    * @function monomerTypeList
    */
    monomerTypeList: function () {
        var monomertypes = { "": "" };
        monomertypes[org.helm.webeditor.HELM.BASE] = "Base";
        monomertypes[org.helm.webeditor.HELM.SUGAR] = "Sugar";
        monomertypes[org.helm.webeditor.HELM.LINKER] = "Linker";
        monomertypes[org.helm.webeditor.HELM.AA] = "Amino Acid";
        monomertypes[org.helm.webeditor.HELM.CHEM] = "Chem";
        return monomertypes;
    },

    symbolCase: function (s) {
        return s == null || org.helm.webeditor.kCaseSensitive ? s : s.toLowerCase();
    },

    /**
    * Show about box
    * @function about
    */
    about: function () {
        var me = this;
        if (this.aboutDlg == null) {
            var div = scil.Utils.createElement(null, "div");
            scil.Utils.createElement(div, "img", null, { width: 425, height: 145 }, { src: scil.Utils.imgSrc("img/helm.png") });

            scil.Utils.createElement(div, "div", "Built on <a target=_blank href='http://www.jsdraw.com'>JSDraw.Lite " + JSDraw2.kFileVersion + "</a> (open source), by <a target=_blank href='http://www.scillignece.com'>Scilligence</a>", { textAlign: "right", paddingRight: "26px" });
            var tbody = scil.Utils.createTable(div, null, null, { borderTop: "solid 1px gray", width: "100%" });
            var tr = scil.Utils.createElement(tbody, "tr");
            scil.Utils.createElement(tr, "td", this.kVersion);
            scil.Utils.createElement(tr, "td", "&copy; 2016, <a target='_blank' href='http://www.pistoiaalliance.org/'>http://www.pistoiaalliance.org/</a>", { textAlign: "center" });
            scil.Utils.createElement(scil.Utils.createElement(tbody, "tr"), "td", "&nbsp;");
            var btn = scil.Utils.createElement(scil.Utils.createElement(div, "div", null, { textAlign: "center" }), "button", "OK", { width: scil.Utils.buttonWidth + "px" });

            me.aboutDlg = new JSDraw2.Dialog("About HELM Web Editor", div);
            scil.connect(btn, "onclick", function (e) { me.aboutDlg.hide(); e.preventDefault(); });
        }
        this.aboutDlg.show();
    },

    isAmbiguous: function (elem, biotype) {
        if (elem == "*" || elem == "_" || biotype == org.helm.webeditor.HELM.AA && elem == 'X' ||
            (biotype == org.helm.webeditor.HELM.SUGAR || biotype == org.helm.webeditor.HELM.BASE || biotype == org.helm.webeditor.HELM.LINKER) && elem == 'N') {
            return true;
        }

        if (!scil.Utils.startswith(elem, '(') || !scil.Utils.endswith(elem, ')'))
            return false;

        elem = elem.substr(1, elem.length - 2);
        var ss = org.helm.webeditor.IO.split(elem, ',');
        if (ss.length > 1)
            return true;

        ss = org.helm.webeditor.IO.split(elem, '+');
        if (ss.length > 1)
            return true;

        return false;
    }
};

scil.helm = org.helm.webeditor;
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* Interface class
* @class org.helm.webeditor.Interface
*/
org.helm.webeditor.Interface = {
    /**
    * Create the canvas
    * @function createCanvas
    * @param {DOM} div
    * @param {dict} args - check <a href='http://www.scilligence.com/sdk/jsdraw/logical/scilligence/JSDraw2/Editor.html'>JSDraw SDK</a>
    */
    createCanvas: function (div, args) {
        return new JSDraw2.Editor(div, args);
    },

    /**
    * Create a molecule object
    * @function createMol
    * @param {string} molfile
    */
    createMol: function (molfile) {
        var m = new JSDraw2.Mol();
        m.setMolfile(molfile);
        return m;
    },

    /**
    * Create a point
    * @function createPoint
    * @param {number} x
    * @param {number} y
    */
    createPoint: function (x, y) {
        return new JSDraw2.Point(x, y);
    },

    /**
    * Create a rectangle object
    * @function createRect
    * @param {number} l - left
    * @param {number} t - top
    * @param {number} w - width
    * @param {number} h - height
    */
    createRect: function (l, t, w, h) {
        return new JSDraw2.Rect(l, t, w, h);
    },

    /**
    * Create an atom
    * @function createAtom
    * @param {JSDraw2.Mol} m
    * @param {JSDraw2.Point} p - the coordinate
    */
    createAtom: function (m, p) {
        return m.addAtom(new JSDraw2.Atom(p));
    },

    /**
    * Create a bond between two atoms
    * @function createBond
    * @param {JSDraw2.Mol} m
    * @param {JSDraw2.Atom} a1
    * @param {JSDraw2.Atom} a2
    */
    createBond: function (m, a1, a2, bondtype) {
        return m.addBond(new JSDraw2.Bond(a1, a2, bondtype == null ? JSDraw2.BONDTYPES.SINGLE : bondtype));
    },

    /**
    * Get atom counts
    * @function getAtomStats
    * @param {JSDraw2.Mol} m
    * @param {array} atoms
    */
    getAtomStats: function (m, atoms) {
        var mol = { atoms: atoms, bonds: m.bonds };
        var ret = JSDraw2.FormulaParser.getAtomStats(m);
        return ret == null ? null : ret.elements;
    },

    /**
    * Test if two molecules are equal
    * @function molEquals
    * @param {JSDraw2.Mol} m1
    * @param {JSDraw2.Mol} m2
    */
    molEquals: function (m1, m2) {
        var mol1 = m1.mol != null ? m1.mol : (m1.mol = this.createMol(scil.helm.Monomers.getMolfile(m1)));
        var mol2 = m2.mol != null ? m2.mol : (m2.mol = this.createMol(scil.helm.Monomers.getMolfile(m2)));
        return mol2.fullstructureMatch(mol1);
    },

    /**
    * count atoms and bonds to calculate MF and MW
    * @function molStats
    * @param {string} molfile
    */
    molStats: function (molfile) {
        var mol = this.createMol(molfile);
        mol.calcHCount();
        return JSDraw2.FormulaParser.getAtomStats(mol).elements;
    },

    /**
    * Get element mass
    * @function getElementMass
    * @param {string} e - element name
    */
    getElementMass: function (e) {
        return JSDraw2.PT[e].m;
    },

    /**
    * Get the current object
    * @function getCurrentAtom
    * @param {JSDraw2.Editor} jsd - JSDraw Editor
    */
    getCurrentAtom: function (jsd) {
        return JSDraw2.Atom.cast(jsd.curObject)
    },

    /**
    * Scale the canvas
    * @function scaleCanvas
    * @param {JSDraw2.Editor} jsd - JSDraw Editor
    */
    scaleCanvas: function (jsd) {
        var scale = JSDraw2.Editor.BONDLENGTH / jsd.bondlength;
        if (JSDraw2.Editor.BONDLENGTH / jsd.bondlength > 1)
            jsd.scale(JSDraw2.Editor.BONDLENGTH / jsd.bondlength);
    },

    /**
    * called by the canvas to draw a monomer
    * @function drawMonomer
    * @param {SVG} surface
    * @param {JSDraw2.Atom} a - monomer object
    * @param {JSDraw2.Point} p - coordinate
    * @param {number} fontsize
    * @param {number} linewidth
    * @param {string} color
    */
    drawMonomer: function (surface, a, p, fontsize, linewidth, color) {
        if (a.hidden)
            return;

        color = null;
        var biotype = a.biotype();
        var c = scil.Utils.isNullOrEmpty(color) ? org.helm.webeditor.Monomers.getColor(a) : color;
        var w = fontsize * org.helm.webeditor.atomscale;
        var lw = linewidth / 2; //(c.nature ? 1 : 2);
        if (biotype == org.helm.webeditor.HELM.LINKER)
            JSDraw2.Drawer.drawEllipse(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.SUGAR)
            JSDraw2.Drawer.drawRect(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw, linewidth * 3).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.BASE)
            JSDraw2.Drawer.drawDiamond(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.AA)
            JSDraw2.Drawer.drawHexgon(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw, linewidth * 3).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.CHEM)
            JSDraw2.Drawer.drawRect(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.BLOB)
            JSDraw2.Drawer.drawRect(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw * 2, linewidth * 5).setFill(c.backgroundcolor);
        else if (biotype == org.helm.webeditor.HELM.NUCLEOTIDE)
            JSDraw2.Drawer.drawPentagon(surface, org.helm.webeditor.Interface.createRect(p.x - w / 2, p.y - w / 2, w, w), c.linecolor, lw, linewidth * 3).setFill(c.backgroundcolor);
        var pt = p.clone();
        p.offset(0, -1);
        JSDraw2.Drawer.drawLabel(surface, p, a.elem, c.textcolor, fontsize * (a.elem.length > 1 ? 2 / a.elem.length : 1.0), null, null, null, false);

        if (a.bio.id > 0) {
            var p1 = p.clone();
            p1.offset(-fontsize * 1.2, -fontsize * 1.2);
            JSDraw2.Drawer.drawLabel(surface, p1, a.bio.id, "#00FF00", fontsize, null, null, null, false);
        }
        if (!scil.Utils.isNullOrEmpty(a.bio.annotation)) {
            var p1 = p.clone();
            var s = a.bio.annotation;
            if (a.bio.annotationshowright) {
                var c = a.biotype() == org.helm.webeditor.HELM.AA ? 0.7 : 1;
                p1.offset(fontsize * c, -fontsize * 1.5);
                JSDraw2.Drawer.drawLabel(surface, p1, s, "#FFA500", fontsize, null, "start", null, false);
            }
            else {
                var c = a.biotype() == org.helm.webeditor.HELM.AA ? 1.5 : 1;
                p1.offset(-fontsize * c, -fontsize * 1.5);
                JSDraw2.Drawer.drawLabel(surface, p1, s, "#FFA500", fontsize, null, "end", null, false);
            }
        }

        if (!scil.Utils.isNullOrEmpty(a.tag)) {
            var r = fontsize / 2;
            JSDraw2.Drawer.drawEllipse(surface, org.helm.webeditor.Interface.createRect(pt.x + r / 2, pt.y - w / 3 - r / 2, r, r), "white", lw).setFill("red");
        }
    },

    addToolbar: function (buttons, flat, sub, options) {
        var sub = [
                { c: "helm_base", t: "Base", label: "Base" },
                { c: "helm_sugar", t: "Sugar", label: "Sugar" },
                { c: "helm_linker", t: "Linker", label: "Linker" },
                { c: "helm_aa", t: "Peptide", label: "Peptide" },
                { c: "helm_chem", t: "Chemistry", label: "Chemistry" }
        ];

        var main = { c: "helm_nucleotide", t: "Nucleotide", label: "Nucleotide", sub: sub, hidden: true };
        buttons.push(main);

        buttons.push({ c: "new", t: "New", label: "New" });
        if (typeof (JSDrawServices) != "undefined" && JSDrawServices.url != null) {
            buttons.push({ c: "open", t: "Load", label: "Load" });
            buttons.push({ c: "save", t: "Save", label: "Save" });
        }
        buttons.push({ c: "|" });
    },

    /**
    * called when the canvas is creating toolbar
    * @function getHelmToolbar
    * @param {array} buttons
    * @param {array} filesubmenus
    * @param {array} selecttools
    * @param {dict} options
    */
    getHelmToolbar: function (buttons, filesubmenus, selecttools, options) {
        this.addToolbar(buttons, true, null, options);

        if (org.helm.webeditor.ambiguity) {
            buttons.push({ c: "helm_blob", t: "BLOB", label: "BLOB" });
            buttons.push({ c: "bracket", t: "Bracket", label: "Bracket" });
        }
        buttons.push({ c: "single", t: "Single bond", label: "Single" });
        buttons.push({ c: "|" });

        buttons.push({ c: "undo", t: "Undo", label: "Undo" });
        buttons.push({ c: "redo", t: "Redo", label: "Redo" });
        buttons.push({ c: "|" });

        buttons.push({ c: "eraser", t: "Eraser", label: "Eraser" });
        buttons.push({ c: "|" });
        buttons.push({ c: "select", t: "Box Selection", label: "Select", sub: selecttools });
        buttons.push({ c: "|" });
        buttons.push({ c: "helm_find", t: "Find/Replace", label: "Find/Replace" });
        buttons.push({ c: "helm_layout", t: "Clean", label: "Clean" });
        buttons.push({ c: "|" });
        buttons.push({ c: "zoomin", t: "Zoom in", label: "Zoom" });
        buttons.push({ c: "zoomout", t: "Zoom out", label: "Zoom" });
        buttons.push({ c: "|" });
        buttons.push({ c: "center", t: "Move to center", label: "Center" });
        buttons.push({ c: "moveview", t: "Move/View", label: "Move" });
    },

    /**
    * called when the canvas is trying to display context menu
    * @function onContextMenu
    * @param {JSDraw2.Editor} ed - JSDraw Editor
    * @param {Event} e - Javascript event
    * @param {bool} viewonly - indicate if this is viewonly mode
    */
    onContextMenu: function (ed, e, viewonly) {
        var items = [];

        if (ed.options.toolbarmode == "helm") {
            var a = JSDraw2.Atom.cast(ed.curObject);
            var b = JSDraw2.Bond.cast(ed.curObject);
            var grp = JSDraw2.Group.cast(ed.curObject);
            if (a != null && a.hidden)
                a = null;
            if (a != null) {
                var biotype = a.biotype();
                if (biotype == scil.helm.HELM.SUGAR && a.bio != null) {
                    items.push({ caption: "Set as Sense", key: "helm_set_sense" });
                    items.push({ caption: "Set as Antisense", key: "helm_set_antisense" });
                    items.push({ caption: "Clear Annotation", key: "helm_set_clear" });
                    items.push("-");
                    items.push({ caption: "Create Complementary Strand", key: "helm_complementary_strand" });
                }

                if (a.superatom != null) {
                    if (items.length > 0)
                        items.push("-");
                    items.push({ caption: "Expand" });
                }

                if (org.helm.webeditor.ambiguity && org.helm.webeditor.isHelmNode(biotype) && a.superatom == null) {
                    if (items.length > 0)
                        items.push("-");
                    if (biotype == org.helm.webeditor.HELM.BLOB)
                        items.push({ caption: "Blob Type", callback: function (cmd, obj) { ed.helm.setHelmBlobType(obj, cmd); }, children: org.helm.webeditor.blobtypes });
                    else if (a.group == null)
                        items.push({ caption: "Create Group", key: "helm_create_group" });
                    items.push("-");
                    items.push({ caption: "Set Monomer Attributes", key: "helm_atom_prop" });
                }
            }
            else if (b != null && (b.a1.biotype() == org.helm.webeditor.HELM.BLOB || b.a2.biotype() == org.helm.webeditor.HELM.BLOB)) {
                items.push({ caption: "Set Bond Attributes", key: "helm_bond_prop" });
            }
            else if (grp != null) {
                items.push({ caption: "Create Group", key: "helm_create_group" });
                items.push("-");
                items.push({ caption: "Collapse", key: "helm_group_collapse" });
                items.push("-");
                items.push({ caption: "Set Attributes", key: "group_setproperties" });
            }
        }
        else {
            var a = JSDraw2.Atom.cast(ed.curObject);
            if (a != null && a.bio == null)
                items.push({ caption: "R Group", callback: function (cmd, obj) { ed.menuSetAtomType(cmd, obj); }, children: ["R1", "R2", "R3", "R4", "R5"] });

            items.push({ caption: "Copy Molfile", key: "copymolfile" });
        }

        var br = JSDraw2.Bracket.cast(ed.curObject);
        if (br != null)
            items.push({ caption: "Set Subscript", key: "setbracketsubscription" });

        if (items.length > 0)
            items.push("-");

        if (ed.options.toolbarmode == "helm")
            ; //items.push({ caption: "About HELM Web Editor", key: "abouthelm" });
        else
            items.push({ caption: "About JSDraw", key: "about" });
        return items;
    }
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* MonomerColors class
* @class org.helm.webeditor.MonomerColors
*/
org.helm.webeditor.MonomerColors = {
    unknown: "#FFFF00",

    bases: {
        A: "#A0A0FF",
        G: "#FF7070",
        T: "#A0FFA0",
        C: "#FF8C4B",
        U: "#FF8080"
    },

    linkers: {
        P: "#9aa5e1",
        p: "#9aa5e1"
    },

    sugars: {
        R: "#7a85c1",
        r: "#7a85c1"
    },

    aas: {
        A: "#C8C8C8",
        R: "#145AFF",
        N: "#00DCDC",
        D: "#E60A0A",
        C: "#E6E600",
        E: "#00DCDC",
        Q: "#E60A0A",
        G: "#EBEBEB",
        H: "#8282D2",
        I: "#0F820F",
        L: "#0F820F",
        K: "#145AFF",
        M: "#E6E600",
        F: "#3232AA",
        P: "#DC9682",
        S: "#FA9600",
        T: "#FA9600",
        W: "#B45AB4",
        Y: "#3232AA",
        V: "#0F820F"
    },

    chems: {
        R: "#eeeeee",
    },

    blobs: {
        B: "#999999",
        G: "#e2e2e2"
    }
};
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/


/**
* Monomers class
* @class org.helm.webeditor.Monomers
*/
org.helm.webeditor.Monomers = {
    smilesmonomerid: 0,
    smilesmonomers: {},
    aliasset: {},
    defaultmonomers: { HELM_BASE: null, HELM_SUGAR: null, HELM_LINKER: null, HELM_AA: null, HELM_CHEM: null },
    blobs: { blob: { n: 'Blob', id: "Blob", na: 'B', rs: 0, at: {}, m: '' }, group: { n: 'Group', id: "Group", na: 'G', rs: 0, at: {}, m: ''} },

    /**
    * Clear monomer database (internal use)
    * @function clear
    */
    clear: function () {
        this.sugars = {};
        this.linkers = {};
        this.bases = {};
        this.aas = {};
        this.chems = {};
    },

    /**
    * Get the default monomer of a given monomer type (internal use)
    * @function getDefaultMonomer
    */
    getDefaultMonomer: function (monomertype) {
        var r = this.defaultmonomers[monomertype];
        if (r != null)
            return r;

        if (monomertype == org.helm.webeditor.HELM.BASE)
            return this._getFirstKey(org.helm.webeditor.Monomers.bases, "a");
        else if (monomertype == org.helm.webeditor.HELM.SUGAR)
            return this._getFirstKey(org.helm.webeditor.Monomers.sugars, "r");
        else if (monomertype == org.helm.webeditor.HELM.LINKER)
            return this._getFirstKey(org.helm.webeditor.Monomers.linkers, "p");
        else if (monomertype == org.helm.webeditor.HELM.AA)
            return this._getFirstKey(org.helm.webeditor.Monomers.aas, "A");
        else if (monomertype == org.helm.webeditor.HELM.CHEM)
            return this._getFirstKey(org.helm.webeditor.Monomers.chems, "R");
        else if (monomertype == org.helm.webeditor.HELM.BLOB)
            return this._getFirstKey(org.helm.webeditor.Monomers.blobs);
        return "?";
    },

    /**
    * Tool function (internal use)
    * @function _getFirstKey
    */
    _getFirstKey: function (set, key1, key2) {
        if (key1 != null && set[scil.helm.symbolCase(key1)] != null)
            return set[scil.helm.symbolCase(key1)].id;
        if (key2 != null && set[scil.helm.symbolCase(key2)] != null)
            return set[scil.helm.symbolCase(key2)].id;

        for (var k in set)
            return k;
        return "?";
    },

    /**
    * Save monomers as text database (internal use)
    * @function saveTextDB
    */
    saveTextDB: function (url) {
        var cols = ["id", "symbol", "name", "naturalanalog", "molfile", "smiles", "polymertype", "monomertype", "r1", "r2", "r3", "r4", "r5", "author", "createddate"];
        var s = "";
        var n = { n: 0 };
        s += this.saveMonomersAsText(this.aas, "PEPTIDE", "Undefined", cols, n);
        s += this.saveMonomersAsText(this.sugars, "RNA", "Backbone", cols, n);
        s += this.saveMonomersAsText(this.linkers, "RNA", "Backbone", cols, n);
        s += this.saveMonomersAsText(this.bases, "RNA", "Branch", cols, n);
        s += this.saveMonomersAsText(this.chems, "CHEM", "Undefined", cols, n);

        s = n.n + "\n" + s;
        if (url == null)
            return s;

        var args = { client: "jsdraw", wrapper: "none", filename: "monomers.txt", directsave: 1, contents: s };
        scil.Utils.post(url, args, "_blank");
    },

    /**
    * Save monomers into xml string (internal use)
    * @function saveMonomerDB
    */
    saveMonomerDB: function (url) {
        var s = "<MONOMER_DB>\n";
        s += "<PolymerList>\n";

        s += "<Polymer polymerType='PEPTIDE'>\n";
        s += this.saveMonomers(this.aas, "PEPTIDE", "Undefined");
        s += "</Polymer>\n";

        s += "<Polymer polymerType='RNA'>\n";
        s += this.saveMonomers(this.sugars, "RNA", "Backbone");
        s += this.saveMonomers(this.linkers, "RNA", "Backbone");
        s += this.saveMonomers(this.bases, "RNA", "Branch");
        s += "</Polymer>\n";

        s += "<Polymer polymerType='CHEM'>\n";
        s += this.saveMonomers(this.chems, "CHEM", "Undefined");
        s += "</Polymer>\n";

        s += "</PolymerList>\n";
        s += "</MONOMER_DB>";

        if (url == null)
            return s;

        var args = { client: "jsdraw", wrapper: "none", filename: "HELMMonomerDB.xml", directsave: 1, contents: s };
        scil.Utils.post(url, args, "_blank");
    },

    /**
    * Save all monomers into a text file
    * @function saveMonomersAsText
    */
    saveMonomersAsText: function (set, type, mt, cols, n) {
        var ret = "";
        for (var id in set) {
            var s = this.writeOneAsText({ id: ++n.n, symbol: id, monomertype: mt, polymertype: type, name: set[id].n, naturalanalog: set[id].na, m: set[id] }, cols);
            ret += JSDraw2.Base64.encode(s) + "\n";
        }

        return ret;
    },

    /**
    * Save all Monomers into xml 
    * @function saveMonomers
    */
    saveMonomers: function (set, type, mt) {
        var s = "";
        for (var id in set)
            s += this.writeOne({ id: id, mt: mt, type: type, m: set[id] });
        return s;
    },

    /**
    * Load monomer from a web service
    * @function loadFromUrl
    */
    loadFromUrl: function (url, callback) {
        var fn = function (xml) {
            org.helm.webeditor.monomers.loadFromXml(xml);
            if (callback != null)
                callback();
        };
        scil.Utils.download(url, fn);
    },

    /**
    * Load monomer from xml string 
    * @function loadFromXml
    */
    loadFromXml: function (s) {
        var doc = scil.Utils.parseXml(s);
        if (doc == null)
            return false;
        this.loadMonomers(doc);
    },

    /**
    * Load monomer from json array coming from database
    * @function loadDB
    */
    loadDB: function (list, makeMon, clearall) {
        if (clearall != false)
            this.clear();

        if (list.length == null && list.list != null)
            list = list.list;

        for (var i = 0; i < list.length; ++i) {
            var x = list[i];

            var m = null;
            if (makeMon != null) {
                m = makeMon(x);
            }
            else {
                m = { id: x.symbol, n: x.name, na: x.naturalanalog, type: x.polymertype, mt: x.monomertype, m: x.molfile };

                m.at = {};
                var rs = 0;
                for (var r = 1; r <= 5; ++r) {
                    if (x["r" + r]) {
                        m.at["R" + r] = x["r" + r];
                        ++rs;
                    }
                }
                m.rs = rs;
            }

            this.addOneMonomer(m);
        }
    },

    /**
    * Load monomer from XML 
    * @function loadMonomers
    */
    loadMonomers: function (doc, callback) {
        var list = doc.getElementsByTagName("Monomer");
        if (list == null || list.length == 0)
            return false;

        if (callback == null) {
            for (var i = 0; i < list.length; ++i) {
                var m = this.readOne(list[i]);
                if (m != null)
                    this.addOneMonomer(m);
            }
            return true;
        }

        var newmonomers = [];
        var overlapped = [];
        for (var i = 0; i < list.length; ++i) {
            var m = this.readOne(list[i]);
            var old = this.getMonomer(this.helm2Type(m), m.id);
            if (old == null)
                newmonomers.push(m);
            else {
                if (!org.helm.webeditor.Interface.molEquals(old, m))
                    overlapped.push(m);
            }
        }

        var me = this;
        this.renameNextMonomer(newmonomers, overlapped, function () {
            var renamed = [];
            for (var i = 0; i < newmonomers.length; ++i) {
                var m = newmonomers[i];
                me.addOneMonomer(m);
                if (m.oldname != null)
                    renamed.push(m);
            }
            callback(renamed);
        });
    },

    /**
    * Rename a monomer (internal use)
    * @function renameNextMonomer
    */
    renameNextMonomer: function (newmonomers, overlapped, callback) {
        if (overlapped.length == 0) {
            callback();
            return;
        }

        var me = this;
        var m = overlapped[0];

        scil.Utils.prompt2({
            caption: "Duplicate Monomer",
            message: "Monomer name, " + m.id + ", is used. Please enter a new name for it:",
            callback: function (s) {
                if (me.getMonomer(m.type, s) == null) {
                    m.oldname = m.id;
                    m.id = s;
                    newmonomers.push(m);
                    overlapped.splice(0, 1);
                }
                me.renameNextMonomer(newmonomers, overlapped, callback);
            }
        });
    },

    getAliases: function (biotype) {
        return null;
    },

    /**
    * Get the monomer set by its type (internal use)
    * @function getMonomerSet
    */
    getMonomerSet: function (a) {
        if (a == null)
            return null;
        if (a.T == "ATOM")
            a = a.biotype();
        if (a == org.helm.webeditor.HELM.BASE)
            return org.helm.webeditor.monomers.bases;
        else if (a == org.helm.webeditor.HELM.SUGAR)
            return org.helm.webeditor.monomers.sugars;
        else if (a == org.helm.webeditor.HELM.LINKER)
            return org.helm.webeditor.monomers.linkers;
        else if (a == org.helm.webeditor.HELM.AA)
            return org.helm.webeditor.monomers.aas;
        else if (a == org.helm.webeditor.HELM.CHEM)
            return org.helm.webeditor.monomers.chems;
        else if (a == org.helm.webeditor.HELM.BLOB)
            return org.helm.webeditor.monomers.blobs;
        return null;
    },

    /**
    * Get all monomer colors (internal use)
    * @function getMonomerColors
    */
    getMonomerColors: function (a) {
        if (a == null)
            return null;
        if (a.T == "ATOM")
            a = a.biotype();
        if (a == org.helm.webeditor.HELM.BASE)
            return org.helm.webeditor.MonomerColors.bases;
        else if (a == org.helm.webeditor.HELM.SUGAR)
            return org.helm.webeditor.MonomerColors.sugars;
        else if (a == org.helm.webeditor.HELM.LINKER)
            return org.helm.webeditor.MonomerColors.linkers;
        else if (a == org.helm.webeditor.HELM.AA)
            return org.helm.webeditor.MonomerColors.aas;
        else if (a == org.helm.webeditor.HELM.CHEM)
            return org.helm.webeditor.MonomerColors.chems;
        else if (a == org.helm.webeditor.HELM.BLOB)
            return org.helm.webeditor.MonomerColors.blobs;
        return null;
    },

    /**
    * Get monomer list of a type (internal use)
    * @function getMonomerList
    */
    getMonomerList: function (a) {
        var set = this.getMonomerSet(a);
        if (set == null)
            return null;

        var ret = [];
        for (var k in set)
            ret.push(set[k].id);

        return ret;
    },

    /**
    * Get a monomer by an object or its name (internal use)
    * @function getMonomer
    */
    getMonomer: function (a, name) {
        if (a == null && name == null)
            return null;

        var s;
        var biotype;
        if (name == null) {
            biotype = a.biotype();
            s = a.elem;
        }
        else {
            biotype = a;
            s = org.helm.webeditor.IO.trimBracket(name);
        }

        if (s == "?") {
            var m = { id: '?', n: "?", na: '?', rs: 2, at: {}, m: "" }; // https://github.com/PistoiaHELM/HELMWebEditor/issues/213
            if (biotype == org.helm.webeditor.HELM.SUGAR)
                m.at.R3 = "H";
            return m;
        }

        if (biotype == org.helm.webeditor.HELM.LINKER && s == "null")
            return { id: 'null', n: "?", na: '?', rs: 2, at: { R1: 'H', R2: 'H' }, m: "" };

        var set = this.getMonomerSet(biotype);
        if (set == null)
            return null;

        var m = set[scil.helm.symbolCase(s)];
        if (m != null)
            return m;

        var set = this.getAliases(biotype);
        if (set == null)
            return null;

        var m = set[scil.helm.symbolCase(s)];
        if (m != null)
            return m;
  },

    /**
    * Check if the monomer have a R group (internal use)
    * @function hasR
    */
    hasR: function (type, name, r) {
        var m = this.getMonomer(type, name);
        return m != null && m.at != null && m.at[r] != null;
    },

    /**
    * Get monomer color by a monomer object (internal use)
    * @function getColor
    */
    getColor: function (a) {
        var m = this.getMonomer(a, a.elem);
        if (m == null)
            m = {};

        var mc = this.getMonomerColors(a);
        if (mc == null)
            mc = {};
        var color = mc[m.na];

        if (m.backgroundcolor == null && a.elem == "?")
            m.backgroundcolor = org.helm.webeditor.MonomerColors.unknown;

        return {
            linecolor: m.linecolor == null ? "#000" : m.linecolor,
            backgroundcolor: m.backgroundcolor == null ? (color == null ? "#eee" : color) : m.backgroundcolor,
            textcolor: m.textcolor == null ? "#000" : m.textcolor,
            nature: m.nature
        };
    },

    /**
    * Get monomer color by type by name (internal use)
    * @function getColor2
    */
    getColor2: function (type, name) {
        var m = this.getMonomer(type, name);
        if (m == null)
            m = {};

        var mc = this.getMonomerColors(type);
        if (mc == null)
            mc = {};
        var color = mc[m.na];

        return {
            linecolor: m.linecolor == null ? "#000" : m.linecolor,
            backgroundcolor: m.backgroundcolor == null ? (color == null ? "#eee" : color) : m.backgroundcolor,
            textcolor: m.textcolor == null ? "#000" : m.textcolor,
            nature: m.nature
        };
    },

    /**
    * Get the molfile of a monomer (internal use)
    * @function getMolfile
    */
    getMolfile: function (m) {
        if (m != null && m.m == null && m.mz != null)
            m.m = org.helm.webeditor.IO.uncompressGz(m.mz);
        return m == null ? null : m.m;
    },

    /**
    * Convert XML type to HELM Editor type (internal use)
    * @function helm2Type
    */
    helm2Type: function (m) {
        if (m.type == "PEPTIDE")
            return org.helm.webeditor.HELM.AA;
        else if (m.type == "CHEM")
            return org.helm.webeditor.HELM.CHEM;
        else if (m.type == "RNA") {
            if (m.mt == "Branch")
                return org.helm.webeditor.HELM.BASE;
            if (m.mt == "Backbone") {
                if (m.na == "P" || m.na == "p")
                    return org.helm.webeditor.HELM.LINKER;
                else
                    return org.helm.webeditor.HELM.SUGAR;
            }
        }
        return null;
    },

    convertAtomMapping2Rs: function (smiles) {
        if (smiles == null)
            return null;

        // represented with atom mapping: CHEM1{[[*:2]C(=O)[C@H](C)N([*:1])CC]}$$$$
        for (var i = 1; i <= 10; ++i)
            smiles = smiles.replace(new RegExp("\\[\\*\\:" + i + "\\]"), "[R" + i + "]");
        return smiles;
    },

    addSmilesMonomer: function (type, smiles) {
        smiles = this.convertAtomMapping2Rs(smiles);
        var ss = this.findSmilesRs(smiles);
        if (ss == null || ss.length == 0)
            return null;

        smiles = this.chemAxon2JSDrawSmiles(smiles);
        if (this.smilesmonomers[smiles] != null)
            return this.smilesmonomers[smiles];

        var m = { at: {}, smiles: smiles, issmiles: true };
        m.id = "#" + (++this.smilesmonomerid);
        m.name = "SMILES Monomer #" + this.smilesmonomerid;
        for (var i = 0; i < ss.length; ++i)
            m.at[ss[i]] = "H";
        m.rs = ss.length;
        var set = this.getMonomerSet(type);
        set[scil.helm.symbolCase(m.id)] = m;

        if (this.cleanupurl != null) {
            if (this.onMonomerSmiles != null) {
                this.onMonomerSmiles(m, smiles);
            }
            else {
                scil.Utils.ajax(this.cleanupurl, function (ret) {
                    if (ret != null && ret.output != null)
                        m.m = ret.output;
                }, { input: smiles, inputformat: "smiles", outputformat: "mol" });
            }
        }

        this.smilesmonomers[smiles] = m;
        return m;
    },

    chemAxon2JSDrawSmiles: function (smiles) {
        return smiles;
    },

    findSmilesRs: function (s) {
        // "C[13C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|"

        var ret = [];
        // JSDraw like Rs
        for (var i = 1; i <= 10; ++i) {
            var s2 = s.replace(new RegExp("\\[R" + i + "\\]"), "");
            if (s2.length == s.length)
                continue;
            s = s2;
            ret.push("R" + i);
        }

        if (ret.length == 0) {
            // ChemAxon like Rs
            for (var i = 1; i <= 10; ++i) {
                var s2 = s.replace(new RegExp("_R" + i), "");
                if (s2.length == s.length)
                    continue;
                s = s2;
                ret.push("R" + i);
            }
        }

        return ret;
    },

    /**
    * add one monomer to HELM Editor (internal use)
    * @function addOneMonomer
    */
    addOneMonomer: function (m) {
        var set = this.getMonomerSet(this.helm2Type(m));
        if (set == null)
            return false;

        delete m.type;
        delete m.mt;

        set[scil.helm.symbolCase(m.id)] = m;
        return true;
    },

    /**
    * Write one monomer into text file (internal use)
    * @function writeOneAsText
    */
    writeOneAsText: function (m, cols) {
        var molfile = m.m.mz;
        if (scil.Utils.isNullOrEmpty(molfile) && m.m.m != null)
            molfile = m.m.m;

        m.molfile = molfile;
        if (m.m.at != null) {
            for (var x in m.m.at)
                m[scil.helm.symbolCase(x)] = m.m.at[x];
        }

        var s = "";
        for (var i = 0; i < cols.length; ++i) {
            if (i > 0)
                s += "|";
            var k = cols[i];
            s += m[k] == null ? "" : m[k];
        }
        return s;
    },

    /**
    * Save one monomer into xml (internal use)
    * @function writeOne
    */
    writeOne: function (m) {
        var molfile = this.getMolfile(m.m);
        if (molfile != null) {
            var s = org.helm.webeditor.IO.compressGz(molfile); // compress molfile
            if (s != null)
                molfile = s;
        }

        var s = "<Monomer>\n";
        s += "<MonomerID>" + scil.Utils.escXmlValue(m.id) + "</MonomerID>\n";
        s += "<MonomerSmiles>" + scil.Utils.escXmlValue(m.smiles) + "</MonomerSmiles>\n";
        s += "<MonomerMolFile>" + scil.Utils.escXmlValue(molfile) + "</MonomerMolFile>\n";
        s += "<NaturalAnalog>" + scil.Utils.escXmlValue(m.m.na) + "</NaturalAnalog>\n";
        s += "<MonomerType>" + scil.Utils.escXmlValue(m.mt) + "</MonomerType>\n";
        s += "<PolymerType>" + scil.Utils.escXmlValue(m.type) + "</PolymerType>\n";
        if (m.m.at != null) {
            s += "<Attachments>\n";
            for (var r in m.m.at) {
                var cap = m.m.at[r];
                s += "<Attachment>\n";
                s += "<AttachmentID>" + r + "-" + cap + "</AttachmentID>\n";
                s += "<AttachmentLabel>" + r + "</AttachmentLabel>\n";
                s += "<CapGroupName>" + cap + "</CapGroupName>\n";
                s += "<CapGroupSmiles></CapGroupSmiles>\n";
                s += "</Attachment>\n";
            }
            s += "</Attachments>\n";
        }
        s += "</Monomer>\n";
        return s;
    },

    /**
    * Read one monomer from XML (internal use)
    * @function readOne
    */
    readOne: function (e) {
        var s = this.readValue(e, "MonomerMolFile");
        var m = null;
        var mz = null;
        if (s != null) {
            if (s.indexOf("M  END") > 0)
                m = s; // uncompressed molfile
            else
                mz = s; // compressed molfile
        }

        var m = {
            type: this.readValue(e, "PolymerType"),
            mt: this.readValue(e, "MonomerType"),
            id: this.readValue(e, "MonomerID"),
            n: this.readValue(e, "MonomerName"),
            na: this.readValue(e, "NaturalAnalog"),
            mz: mz,
            m: m,
            at: {}
        };

        var rs = 0;
        var list = e.getElementsByTagName("Attachment");
        if (list != null) {
            for (var i = 0; i < list.length; ++i) {
                var a = list[i];
                var r = this.readValue(a, "AttachmentLabel");
                var cap = this.readValue(a, "CapGroupName");
                if (m.at[r] == null)
                    ++rs;
                m.at[r] = cap;
            }
        }

        m.rs = rs;
        return m;
    },

    /**
    * Tool function to ready XML text (internal use)
    * @function readValue
    */
    readValue: function (e, name) {
        var list = e.getElementsByTagName(name);
        if (list == null || list.length == 0)
            return null;
        return scil.Utils.getInnerText(list[0]);
    }
};

org.helm.webeditor.monomers = org.helm.webeditor.Monomers;



scil.helm.Monomers.aas = {
    'd': { id: 'D', n: 'Aspartic acid', na: 'D', m: '\n  Marvin  12021015502D          \n\n 10  9  0  0  0  0            999 V2000\n    7.4800   -2.7926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.7655   -3.2050    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    6.0510   -2.7925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.3365   -3.2049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6220   -2.7924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.3364   -4.0299    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    7.4801   -1.9676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.7654   -4.0300    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1945   -3.2052    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    6.0510   -4.4425    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  7  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  8  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4  5  2  0  0  0  0\n  8 10  1  0  0  0  0\nM  RGP  3   6   3   9   2  10   1\nM  END\n\n$$$$\n', rs: 3, at: { R2: 'OH', R3: 'OH', R1: 'H'} },
    'e': { id: 'E', n: 'Glutamic acid', na: 'E', m: '\n  Marvin  12021015512D          \n\n 11 10  0  0  0  0            999 V2000\n    0.8787   -0.2594    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    0.1643    0.1532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5503   -0.2592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2646    0.1534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5933    0.1530    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9792   -0.2590    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2645    0.9784    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5934    0.9780    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8786   -1.0844    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3077   -0.2596    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1641   -1.4969    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  5  1  1  0  0  0  0\n  1  9  1  0  0  0  0\n  1  2  1  1  0  0  0\n  2  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  7  2  0  0  0  0\n  4  6  1  0  0  0  0\n  5  8  2  0  0  0  0\n  5 10  1  0  0  0  0\n  9 11  1  0  0  0  0\nM  RGP  3   6   3  10   2  11   1\nM  END\n\n$$$$\n', rs: 3, at: { R2: 'OH', R3: 'OH', R1: 'H'} },
    'f': { id: 'F', n: 'Phenylalanine', na: 'F', m: '\n  Marvin  08190815502D          \n\n 13 13  0  0  0  0            999 V2000\n   -3.6075    2.0774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.3219    1.6650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.3220    0.8400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.6077    0.4274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8931    0.8398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8930    1.6648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1785    2.0772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4640    1.6646    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -1.4641    0.8396    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7495    2.0770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7493    2.9021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0351    1.6644    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1786    0.4271    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  6  2  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  2  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  0  0  0  0\n  8  7  1  1  0  0  0\n  8  9  1  0  0  0  0\n  8 10  1  0  0  0  0\n 10 11  2  0  0  0  0\n 10 12  1  0  0  0  0\n  9 13  1  0  0  0  0\nM  RGP  2  12   2  13   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'g': { id: 'G', n: 'Glycine', na: 'G', m: '\n  Marvin  08190815292D          \n\n  6  5  0  0  0  0            999 V2000\n   -0.7189   -0.6069    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0044   -0.1945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0043    0.6304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7190   -1.4318    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7101   -0.6071    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4335   -1.8443    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  4  1  1  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  2  5  1  0  0  0  0\n  4  6  1  0  0  0  0\nM  RGP  2   5   2   6   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'a': { id: 'A', n: 'Alanine', na: 'A', m: '\n  Marvin  06250814262D          \n\n  7  6  0  0  0  0            999 V2000\n    5.4886   -3.0482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.2031   -3.4608    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    6.9176   -3.0483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.2030   -4.2858    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    6.9177   -2.2233    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    7.6321   -3.4609    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    5.4886   -4.6983    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1  0  0  0\n  4  2  1  0  0  0  0\n  2  3  1  0  0  0  0\n  3  5  2  0  0  0  0\n  3  6  1  0  0  0  0\n  4  7  1  0  0  0  0\nM  RGP  2   6   2   7   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'c': { id: 'C', n: 'Cysteine', na: 'C', m: '\n  Marvin  12021015502D          \n\n  9  8  0  0  0  0            999 V2000\n    0.8481    0.4124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1336    0.0000    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -0.5808    0.4126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2952    0.0001    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1335   -0.8249    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8482    1.2374    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5625   -0.0001    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5810   -1.2374    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0097    0.4126    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  6  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1  7  1  0  0  0  0\n  2  5  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  5  8  1  0  0  0  0\n  4  9  1  0  0  0  0\nM  RGP  3   7   2   8   1   9   3\nM  END\n\n$$$$\n', rs: 3, at: { R2: 'OH', R1: 'H', R3: 'H'} },
    'l': { id: 'L', n: 'Leucine', na: 'L', m: '\n  Marvin  08200815002D          \n\n 10  9  0  0  0  0            999 V2000\n   -2.7541    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4686    0.9102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4686    1.7352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.1830    2.1477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0397    1.7351    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -1.3250    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3250    2.9726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0397    0.9101    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6105    1.7350    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7542    0.4976    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  3  1  1  0  0  0  0\n  5  1  1  1  0  0  0\n  3  2  1  0  0  0  0\n  3  4  1  0  0  0  0\n  5  8  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  2  0  0  0  0\n  6  9  1  0  0  0  0\n  8 10  1  0  0  0  0\nM  RGP  2   9   2  10   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'm': { id: 'M', n: 'Methionine', na: 'M', m: '\n  Marvin  08190815482D          \n\n 10  9  0  0  0  0            999 V2000\n   -4.1128    0.8183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.8273   -0.4191    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.1127    1.6433    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.8272    0.4059    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -5.5418    0.8185    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2563    0.4061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.6852    0.4062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.9707    0.8186    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3983    0.4057    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -5.5418   -0.8316    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  3  1  2  0  0  0  0\n  1  4  1  0  0  0  0\n  1  9  1  0  0  0  0\n  4  2  1  0  0  0  0\n  4  5  1  1  0  0  0\n  5  6  1  0  0  0  0\n  6  8  1  0  0  0  0\n  8  7  1  0  0  0  0\n  2 10  1  0  0  0  0\nM  RGP  2   9   2  10   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'n': { id: 'N', n: 'Asparagine', na: 'N', m: '\n  Marvin  08190815142D          \n\n 10  9  0  0  0  0            999 V2000\n    6.9268   -2.8241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.2123   -3.2366    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    5.4978   -2.8241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.9268   -1.9991    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.2123   -4.0616    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7834   -3.2366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0690   -2.8242    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7834   -4.0616    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    7.6413   -3.2366    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    5.6290   -4.6450    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  4  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  5  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  6  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6  7  2  0  0  0  0\n  5 10  1  0  0  0  0\nM  RGP  2   9   2  10   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'h': { id: 'H', n: 'Histidine', na: 'H', m: '\n  Marvin  08190815312D          \n\n 12 12  0  0  0  0            999 V2000\n    1.2557    0.4210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5412    0.0086    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -0.1733    0.4211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8878    0.0087    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1685   -0.3674    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9114   -0.8159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6648    0.2860    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7029   -1.0484    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5411   -0.8164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2558    1.2461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9702    0.0084    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1733   -1.2289    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n 10  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1 11  1  0  0  0  0\n  2  9  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  6  4  2  0  0  0  0\n  7  4  1  0  0  0  0\n  7  5  2  0  0  0  0\n  5  8  1  0  0  0  0\n  8  6  1  0  0  0  0\n  9 12  1  0  0  0  0\nM  RGP  2  11   2  12   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'i': { id: 'I', n: 'Isoleucine', na: 'I', m: '\n  Marvin  08190815422D          \n\n 10  9  0  0  0  0            999 V2000\n   -0.7783   -0.6153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4928   -0.2028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6506   -0.6155    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    1.3652   -0.2031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3653    0.6219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6505   -1.4405    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0638   -0.2029    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -0.0637    0.6221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0796   -0.6157    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0640   -1.8530    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  1  0  0  0  0\n  1  2  1  0  0  0  0\n  3  7  1  0  0  0  0\n  3  6  1  0  0  0  0\n  3  4  1  1  0  0  0\n  4  5  2  0  0  0  0\n  4  9  1  0  0  0  0\n  7  8  1  6  0  0  0\n  6 10  1  0  0  0  0\nM  RGP  2   9   2  10   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'k': { id: 'K', n: 'Lysine', na: 'K', m: '\n  Marvin  12021015522D          \n\n 12 11  0  0  0  0            999 V2000\n   -1.2102    2.0919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9248    1.6794    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -2.6392    2.0920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3538    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0682    2.0920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7828    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.4972    2.0921    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9248    0.8544    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2102    2.9169    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4958    1.6794    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6393    0.4419    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2117    1.6796    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  9  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1 10  1  0  0  0  0\n  2  8  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  0  0  0  0\n  8 11  1  0  0  0  0\n  7 12  1  0  0  0  0\nM  RGP  3  10   2  11   1  12   3\nM  END\n\n$$$$\n', rs: 3, at: { R2: 'OH', R1: 'H', R3: 'H'} },
    't': { id: 'T', n: 'Threonine', na: 'T', m: '\n  Marvin  08190820372D          \n\n  9  8  0  0  0  0            999 V2000\n   -3.2927    2.1067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0072    0.8691    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.2926    2.9316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0071    1.6941    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -5.4360    1.6943    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7216    2.1068    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.7215    2.9317    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.5783    1.6942    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7217    0.4566    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  3  1  2  0  0  0  0\n  4  1  1  1  0  0  0\n  1  8  1  0  0  0  0\n  4  2  1  0  0  0  0\n  4  6  1  0  0  0  0\n  6  5  1  0  0  0  0\n  6  7  1  1  0  0  0\n  2  9  1  0  0  0  0\nM  RGP  2   8   2   9   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'w': { id: 'W', n: 'Tryptophan', na: 'W', m: '\n  Marvin  08190820412D          \n\n 16 17  0  0  0  0            999 V2000\n   -0.4698    2.4303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4698    3.2553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1843    1.1927    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1843    2.0177    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -1.8988    2.4303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6133    2.0177    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6271    1.1929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3935    2.2860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.8895    1.6267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.7164    3.0451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5353    3.1451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0313    2.4859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7084    1.7267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4159    0.9513    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2447    2.0178    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8988    0.7802    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  2  0  0  0  0\n  1  4  1  0  0  0  0\n  1 15  1  0  0  0  0\n  4  3  1  0  0  0  0\n  4  5  1  1  0  0  0\n  5  6  1  0  0  0  0\n  7  6  2  0  0  0  0\n  8  6  1  0  0  0  0\n 14  7  1  0  0  0  0\n  8  9  2  0  0  0  0\n  8 10  1  0  0  0  0\n  9 14  1  0  0  0  0\n 13  9  1  0  0  0  0\n 10 11  2  0  0  0  0\n 11 12  1  0  0  0  0\n 12 13  2  0  0  0  0\n  3 16  1  0  0  0  0\nM  RGP  2  15   2  16   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'v': { id: 'V', n: 'Valine', na: 'V', m: '\n  Marvin  08190820502D          \n\n  9  8  0  0  0  0            999 V2000\n   -5.2963    0.6040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.0107    0.1915    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -6.7252    1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.7252    0.6040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.4396    0.1915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.2963    1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.0107   -0.6335    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5818    0.1915    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -6.7252   -1.0460    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  6  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1  8  1  0  0  0  0\n  2  7  1  0  0  0  0\n  2  4  1  1  0  0  0\n  4  3  1  0  0  0  0\n  4  5  1  0  0  0  0\n  7  9  1  0  0  0  0\nM  RGP  2   8   2   9   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'q': { id: 'Q', n: 'Glutamine', na: 'Q', m: '\n  Marvin  08190815272D          \n\n 11 10  0  0  0  0            999 V2000\n    2.8012   -4.0307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0867   -4.4433    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    1.3722   -4.0307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6578   -4.4433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8013   -3.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7711   -4.4432    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0867   -5.2683    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0565   -3.2056    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0566   -4.0306    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5157   -4.4432    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    1.3722   -5.6808    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  5  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1 10  1  0  0  0  0\n  2  7  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  4  9  1  0  0  0  0\n  9  6  1  0  0  0  0\n  9  8  2  0  0  0  0\n  7 11  1  0  0  0  0\nM  RGP  2  10   2  11   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'p': { id: 'P', n: 'Proline', na: 'P', m: '\n  Marvin  08190815522D          \n\n  9  9  0  0  0  0            999 V2000\n   -3.4704    1.8073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9525    1.1377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4647    0.4725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6812    0.7308    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6848    1.5558    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -1.9703    1.9683    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9703    2.7933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2558    1.5558    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0978    0.1474    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  5  1  1  1  0  0  0\n  2  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  2  0  0  0  0\n  6  8  1  0  0  0  0\n  4  9  1  0  0  0  0\nM  RGP  2   8   2   9   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    's': { id: 'S', n: 'Serine', na: 'S', m: '\n  Marvin  08190820582D          \n\n  8  7  0  0  0  0            999 V2000\n   -2.0394    2.5340    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0394    3.3589    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7539    1.2964    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7539    2.1214    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -3.4683    2.5340    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.1827    2.1215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3249    2.1215    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4684    0.8839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  2  0  0  0  0\n  1  4  1  0  0  0  0\n  1  7  1  0  0  0  0\n  4  3  1  0  0  0  0\n  4  5  1  1  0  0  0\n  5  6  1  0  0  0  0\n  3  8  1  0  0  0  0\nM  RGP  2   7   2   8   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'r': { id: 'R', n: 'Arginine', na: 'R', m: '\n  Marvin  08190814412D          \n\n 13 12  0  0  0  0            999 V2000\n    8.4822   -2.8352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7677   -3.2478    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    7.0532   -2.8351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3388   -3.2477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6244   -2.8351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7676   -4.0214    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    8.4822   -2.0102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9099   -3.2477    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1954   -2.8351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1954   -2.0101    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4810   -3.2477    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    9.1967   -3.2477    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    6.9708   -4.2349    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  7  1  2  0  0  0  0\n  1  2  1  0  0  0  0\n  1 12  1  0  0  0  0\n  2  6  1  0  0  0  0\n  2  3  1  1  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  8  1  0  0  0  0\n  8  9  1  0  0  0  0\n  9 11  1  0  0  0  0\n  9 10  2  0  0  0  0\n  6 13  1  0  0  0  0\nM  RGP  2  12   2  13   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'} },
    'y': { id: 'Y', n: 'Tyrosine', na: 'Y', m: '\n  Marvin  08190820432D          \n\n 14 14  0  0  0  0            999 V2000\n   -2.8122    1.4277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5268    0.1903    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8121    2.2527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5267    1.0153    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.2412    1.4279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.9557    1.0154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.6701    1.4280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.3846    1.0155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.6702   -0.2220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.9558    0.1904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.3847    0.1905    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.0992   -0.2219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0977    1.0151    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -4.2412   -0.2222    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  3  1  2  0  0  0  0\n  1  4  1  0  0  0  0\n  1 13  1  0  0  0  0\n  4  2  1  0  0  0  0\n  4  5  1  1  0  0  0\n  5  6  1  0  0  0  0\n 10  6  1  0  0  0  0\n  7  6  2  0  0  0  0\n  7  8  1  0  0  0  0\n  8 11  2  0  0  0  0\n 11  9  1  0  0  0  0\n  9 10  2  0  0  0  0\n 11 12  1  0  0  0  0\n  2 14  1  0  0  0  0\nM  RGP  2  13   2  14   1\nM  END\n\n$$$$\n', rs: 2, at: { R2: 'OH', R1: 'H'}}
};
scil.helm.Monomers.sugars = {
    'r': { id: 'R', n: 'Ribose', na: 'R', m: '\n  Marvin  06150820452D          \n\n 12 12  0  0  0  0            999 V2000\n    1.4617    2.2807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7479    1.8672    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    2.1768    1.8692    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    0.9625    1.0707    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n    1.9644    1.0721    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    0.9637    0.2457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9656    0.2471    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5994    2.7216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1756    2.9004    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1279    3.0752    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6971    3.2298    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    1.6650   -0.3083    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  2  4  1  0  0  0  0\n  2  8  1  1  0  0  0\n  3  5  1  0  0  0  0\n  3  9  1  1  0  0  0\n  4  5  1  0  0  0  0\n  4  6  1  6  0  0  0\n  5  7  1  6  0  0  0\n  6 12  1  0  0  0  0\n  8 10  1  0  0  0  0\n 10 11  1  0  0  0  0\nM  RGP  3   9   3  11   1  12   2\nM  END\n\n$$$$\n', rs: 3, at: { R3: 'OH', R1: 'H', R2: 'H'}}
};
scil.helm.Monomers.linkers = {
    'p': { id: 'P', n: 'Phosphate', na: 'P', m: '\n  Marvin  06150820472D          \n\n  5  4  0  0  0  0            999 V2000\n    0.1179    0.7366    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7071    0.7366    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1179    1.5616    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9429    0.7366    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1179   -0.0884    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  RGP  2   2   1   4   2\nM  END\n\n$$$$\n', rs: 2, at: { R1: 'OH', R2: 'OH'}}
};
scil.helm.Monomers.bases = {
    'g': { id: 'G', n: 'Guanine', na: 'G', m: '\n  Marvin  06150820502D          \n\n 12 13  0  0  0  0            999 V2000\n   -1.7509   -3.8489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0364   -4.2614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0364   -5.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7509   -5.4989    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4654   -5.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4654   -4.2614    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3220   -5.4989    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3925   -5.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3925   -4.2614    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7509   -3.0239    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1798   -5.4989    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3220   -6.3239    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1 10  2  0  0  0  0\n  1  6  1  0  0  0  0\n  1  2  1  0  0  0  0\n  9  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  7  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  2  0  0  0  0\n  5  6  1  0  0  0  0\n  5 11  1  0  0  0  0\n  7  8  1  0  0  0  0\n  7 12  1  0  0  0  0\n  8  9  2  0  0  0  0\nM  RGP  1  12   1\nM  END\n\n$$$$\n', rs: 1, at: { R1: 'H'} },
    'a': { id: 'A', n: 'Adenine', na: 'A', m: '\n  Marvin  06150816552D          \n\n 11 12  0  0  0  0            999 V2000\n   -4.0674   -0.7795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3529   -1.1920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3529   -2.0170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0674   -2.4295    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7818   -2.0170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7818   -1.1920    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6384   -2.4295    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9240   -2.0170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9240   -1.1920    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0674    0.0455    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6384   -3.2545    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1 10  1  0  0  0  0\n  1  6  2  0  0  0  0\n  1  2  1  0  0  0  0\n  9  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  7  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  2  0  0  0  0\n  5  6  1  0  0  0  0\n  7  8  1  0  0  0  0\n  7 11  1  0  0  0  0\n  8  9  2  0  0  0  0\nM  RGP  1  11   1\nM  END\n\n$$$$\n', rs: 1, at: { R1: 'H'} },
    'c': { id: 'C', n: 'Cytosine', na: 'C', m: '\n  Marvin  06150820522D          \n\n  9  9  0  0  0  0            999 V2000\n   -1.3611   -7.5896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6466   -8.0021    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6466   -8.8271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3611   -9.2396    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0755   -8.8271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0755   -8.0021    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3611   -6.7646    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7900   -9.2396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3611  -10.0646    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  6  2  0  0  0  0\n  1  7  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  4  9  1  0  0  0  0\n  5  6  1  0  0  0  0\n  5  8  2  0  0  0  0\nM  RGP  1   9   1\nM  END\n\n$$$$\n', rs: 1, at: { R1: 'H'} },
    'u': { id: 'U', n: 'Uracil', na: 'U', m: '\n  Marvin  06150820512D          \n\n  9  9  0  0  0  0            999 V2000\n   -1.2138  -11.7383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4993  -12.1508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4993  -12.9758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2138  -13.3884    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9283  -12.9758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9283  -12.1508    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2138  -10.9133    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6427  -13.3883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2138  -14.2134    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  2  0  0  0  0\n  1  6  1  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  4  9  1  0  0  0  0\n  5  8  2  0  0  0  0\n  5  6  1  0  0  0  0\nM  RGP  1   9   1\nM  END\n\n$$$$\n', rs: 1, at: { R1: 'H'} },
    't': { id: 'T', n: 'Thymine', na: 'T', m: '\n  Marvin  06150820542D          \n\n 10 10  0  0  0  0            999 V2000\n    5.8512   -2.0140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5686   -2.4215    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5744   -3.2465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.8628   -3.6641    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1454   -3.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1396   -2.4315    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    5.8455   -1.1891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.4339   -3.6740    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.8686   -4.4891    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    7.3655   -2.2080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  2  0  0  0  0\n  1  6  1  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  4  9  1  0  0  0  0\n  5  8  2  0  0  0  0\n  5  6  1  0  0  0  0\n  2 10  1  0  0  0  0\nM  RGP  1   9   1\nM  END\n\n$$$$\n', rs: 1, at: { R1: 'H'}}
};
scil.helm.Monomers.chems = {
    'example': { id: 'Example', n: 'Symmetric Doubler', na: null, m: '\n  Marvin  09241011262D          \n\n 23 22  0  0  0  0            999 V2000\n   -3.8304    2.5045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.8304    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1159    1.2670    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4014    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6869    1.2670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9725    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2580    1.2670    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4565    1.6795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4565    2.5045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6869    0.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1709    2.9170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1709    3.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8854    4.1545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8854    4.9795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5448    2.9170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5448    3.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.2593    4.1545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.2593    4.9795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1709    1.2670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5448    1.2670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0994   -0.2725    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -5.9738    5.3920    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    2.5999    5.3920    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  0  0  0  0\n  7  8  1  0  0  0  0\n  8  9  1  0  0  0  0\n  5 10  1  0  0  0  0\n  9 11  1  0  0  0  0\n 11 12  1  0  0  0  0\n 12 13  1  0  0  0  0\n 13 14  1  0  0  0  0\n  1 15  1  0  0  0  0\n 15 16  1  0  0  0  0\n 16 17  1  0  0  0  0\n 17 18  1  0  0  0  0\n  8 19  2  0  0  0  0\n  2 20  2  0  0  0  0\n 10 21  1  0  0  0  0\n 18 22  1  0  0  0  0\n 14 23  1  0  0  0  0\nM  RGP  3  21   1  22   2  23   3\nM  END\n\n$$$$\n', rs: 3, at: { R1: 'H', R2: 'H', R3: 'H'} },
    'r': { id: 'R', n: 'R', na: null, m: null, rs: 0, at: {}}
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* HELM Editor Plugin class
* @class org.helm.webeditor.Plugin
*/
org.helm.webeditor.Plugin = scil.extend(scil._base, {
    /**
    @property {MonomerExplorer} monomerexplorer - Monomer Explorer
    **/
    /**
    @property {JSDraw2.Editor} jsd - Drawing Canvas
    **/

    /**
    * @constructor Plugin
    * @param {JSDraw2.Editor} jsd - The JSDraw canvas
    **/
    constructor: function (jsd) {
        this.jsd = jsd;
        this.monomerexplorer = null;
    },

    /**
    * Get the molecule formula
    * @function getMF
    * @param {bool} html - indicate if html format is needed
    * @returns the molecular formula as a string
    */
    getMF: function (html) {
        return org.helm.webeditor.Formula.getMF(this.jsd.m, html);
    },

    /**
    * Get the molecule weight
    * @function getMW
    * @returns the molecular weight as a number
    */
    getMW: function () {
        return org.helm.webeditor.Formula.getMW(this.jsd.m);
    },

    /**
    * Get the Extinction Coefficient
    * @function getExtinctionCoefficient
    * @returns the Extinction Coefficient as a number
    */
    getExtinctionCoefficient: function () {
        return org.helm.webeditor.ExtinctionCoefficient.calculate(this.jsd.m);
    },

    getSpareRs: function (a, rs) {
        if (a.bio == null) // not bio
            return [];

        var m = org.helm.webeditor.Monomers.getMonomer(a);
        if (m == null)
            return null;

        // https://github.com/PistoiaHELM/HELMWebEditor/issues/213
        if (m.id == "?")
            return "?";

        if (rs == null)
            rs = [];
        else
            rs.splice(0, rs.length);

        for (var r in m.at) {
            var i = parseInt(r.substr(1));
            rs[i] = true;
        }

        var used = [];
        var bonds = this.jsd.m.getNeighborBonds(a);
        for (var i = 0; i < bonds.length; ++i) {
            var b = bonds[i];
            if (b.a1 == a) {
                used[b.r1] = true;
                if (rs[b.r1] != null)
                    rs[b.r1] = false;
            }
            else if (b.a2 == a) {
                used[b.r2] = true;
                if (rs[b.r2] != null)
                    rs[b.r2] = false;
            }
        }

        var ret = [];
        for (var i = 1; i <= rs.length; ++i) {
            if (rs[i])
                ret.push(i);
        }

        if (ret.length == 0 && a.biotype() == org.helm.webeditor.HELM.BLOB)
            return "?";

        return ret.length == 0 ? null : ret;
    },

    setAtomProp: function (obj) {
        var a = JSDraw2.Atom.cast(obj);
        if (a == null)
            return;

        if (this.atomPropDlg == null) {
            var me = this;
            var fields = {
                elem: { label: "Monomer Symbol" },
                tag: { label: "Annotation" }
            };
            this.atomPropDlg = scil.Form.createDlgForm("Monomer Attributes", fields, { label: "Save", onclick: function () { me.setAtomProp2(); } });
        }
        else {
            this.atomPropDlg.show();
        }

        var elem = a.elem;
        if (elem == "?" && !scil.Utils.isNullOrEmpty(a.bio.ambiguity))
            elem = a.bio.ambiguity;

        if (elem == "Blob" || elem == "Group")
            this.atomPropDlg.form.fields.elem.setAttribute("readonly", "readonly");
        else
            this.atomPropDlg.form.fields.elem.removeAttribute("readonly");

        var data = { elem: elem, tag: a.tag };
        this.atomPropDlg.form.setData(data);
        this.atomPropDlg.atom = a;
    },

    setAtomProp2: function () {
        var data = this.atomPropDlg.form.getData();
        data.elem = scil.Utils.trim(data.elem);
        data.tag = scil.Utils.trim(data.tag);

        if (scil.Utils.isNullOrEmpty(data.elem)) {
            scil.Utils.alert("Monomer Type cannot be blank");
            return;
        }

        var f = false;
        var clone = this.jsd.clone();
        var a = this.atomPropDlg.atom;
        if (a.elem != data.elem || (a.tag == null ? "" : a.tag) != data.tag) {
            var set = org.helm.webeditor.Monomers.getMonomerSet(a.biotype());
            var m = set[scil.helm.symbolCase(data.elem)];
            if (m == null) {
                if (!org.helm.webeditor.isAmbiguous(data.elem, a.biotype())) {
                    data.elem = "(" + data.elem + ")";
                    if (!org.helm.webeditor.isAmbiguous(data.elem, a.biotype())) {
                        scil.Utils.alert("Invalid Monomer Type");
                        return;
                    }
                }

                a.elem = "?";
                a.bio.ambiguity = data.elem;
            }
            else {
                a.elem = data.elem;
                a.bio.ambiguity = null;
            }

            f = true;
            a.tag = data.tag == "" ? null : data.tag;
        }

        this.atomPropDlg.hide();

        if (f) {
            this.jsd.pushundo(clone);
            this.jsd.refresh(true);
        }
    },

    setBondProp: function (obj) {
        var b = JSDraw2.Bond.cast(obj);
        if (b == null)
            return;

        var data = {};
        var blob1 = this._getDisplayName(b, data, 1);
        var blob2 = this._getDisplayName(b, data, 2);
        if (!blob1 && !blob2)
            return;

        data.a1ratio = b.ratio1;
        data.a2ratio = b.ratio2;

        if (this.bondPropDlg == null) {
            var me = this;
            var fields = {
                a1: { label: "Monomer #1", viewonly: true },
                a1pos: { label: "Position" },
                a1r: { label: "R#" },
                a1ratio: { label: "Ratio", type: "number", accepts: "^[?]$" },
                a2: { label: "Monomer #2", viewonly: true },
                a2pos: { label: "Position" },
                a2r: { label: "R#" },
                a2ratio: { label: "Ratio", type: "number", accepts: "^[?]$" }
            };
            this.bondPropDlg = scil.Form.createDlgForm("Bond Attributes", fields, { label: "Save", onclick: function () { me.setBondProp2(); } });
        }
        else {
            this.bondPropDlg.show();
        }

        this.bondPropDlg.form.fields.a1pos.disabled = !blob1;
        this.bondPropDlg.form.fields.a1r.disabled = !blob1;
        this.bondPropDlg.form.fields.a2pos.disabled = !blob2;
        this.bondPropDlg.form.fields.a2r.disabled = !blob2;

        if (scil.Utils.isNullOrEmpty(data.a1ratio))
            data.a1ratio = org.helm.webeditor.defaultbondratio;
        if (scil.Utils.isNullOrEmpty(data.a2ratio))
            data.a2ratio = org.helm.webeditor.defaultbondratio;

        this.bondPropDlg.form.setData(data);
        this.bondPropDlg.bond = b;
    },

    _getDisplayName: function (b, data, i) {
        var ai = "a" + i;
        var ri = "r" + i;
        var air = "a" + i + "r";
        var aipos = "a" + i + "pos";

        var ret = false;
        if (b[ai].biotype() == org.helm.webeditor.HELM.BLOB) {
            ret = true;
            data[ai] = b[ai].elem;
            var p = typeof (b[ri]) == "string" ? b[ri].indexOf(':') : -1;
            if (p > 0) {
                data[aipos] = b[ri].substr(0, p);
                data[air] = b[ri].substr(p + 1);
            }
            else {
                data[air] = b[ri];
            }
        }
        else {
            data[ai] = b[ai].elem;
            data.a1r = b[ri];
        }

        return ret;
    },

    setBondProp2: function () {
        var data = this.bondPropDlg.form.getData();
        this.bondPropDlg.hide();

        if (data.a1ratio == "?" || data.a2ratio == "?") {
            data.a1ratio = data.a2ratio = "?";
        }
        else {
            if (data.a1ratio == "")
                data.a1ratio = org.helm.webeditor.defaultbondratio;
            if (data.a2ratio == "")
                data.a2ratio = org.helm.webeditor.defaultbondratio;
        }

        var f = false;
        var clone = this.jsd.clone();
        var b = this.bondPropDlg.bond;
        if (b.a1.biotype() == org.helm.webeditor.HELM.BLOB) {
            if (this._makeBondR(data.a1pos, data.a1r, b, "r1"))
                f = true;
        }
        if (!this._isRatioEq(b.ratio1, data.a1ratio)) {
            b.ratio1 = data.a1ratio > 0 || data.a1ratio == "?" ? data.a1ratio : null;
            f = true;
        }

        if (b.a2.biotype() == org.helm.webeditor.HELM.BLOB) {
            if (this._makeBondR(data.a2pos, data.a2r, b, "r2"))
                f = true;
        }
        if (!this._isRatioEq(b.ratio2, data.a2ratio)) {
            b.ratio2 = data.a2ratio > 0 || data.a2ratio == "?" ? data.a2ratio : null;
            f = true;
        }

        if (f) {
            this.jsd.pushundo(clone);
            this.jsd.refresh(true);
        }
    },

    _isRatioEq: function (r1, r2) {
        if (!(r1 > 0 || r1 == "?") || r1 == "")
            r1 = null;
        if (!(r2 > 0 || r2 == "?") || r2 == "")
            r2 = null;

        return r1 == r2;
    },

    _makeBondR: function (pos, r, b, key) {
        if (pos != null)
            pos = pos.trim();
        if (scil.Utils.isNullOrEmpty(pos))
            pos = "?";

        if (r != null)
            r = r.trim();
        if (scil.Utils.isNullOrEmpty(r))
            r = "?";

        var s = pos + ":" + (!isNaN(r) && r > 0 ? "R" + r : r);
        if (b[key] == s)
            return false;

        b[key] = s;
        return true;
    },

    expandSuperAtom: function (obj) {
        var a = JSDraw2.Atom.cast(obj);
        if (a == null || a.bio == null || !org.helm.webeditor.isHelmNode(a))
            return false;

        this.jsd.pushundo();
        if (this.groupExpand(a))
            this.jsd.refresh(true);
        return true;
    },

    hasSpareR: function (a, r) {
        if (a == null)
            return false;
        if (a.bio == null)
            return true;

        if (typeof (r) == "string" && scil.Utils.startswith(r, "R"))
            r = parseInt(r.substr(1));

        if (a.biotype() == org.helm.webeditor.HELM.BLOB)
            return true;

        var rs = this.getSpareRs(a);
        if (rs == null || rs.indexOf(r) < 0) {
            //scil.Utils.alert("The monomer, " + a.elem + ", does define R" + r);
            return false;
        }

        var bonds = this.jsd.m.getNeighborBonds(a);
        for (var i = 0; i < bonds.length; ++i) {
            var b = bonds[i];
            if (b.a1 == a && b.r1 == r)
                return false;
            else if (b.a2 == a && b.r2 == r)
                return false;
        }

        return true;
    },

    getDefaultNodeType: function (a, c) {
        var s = null;
        if (this.monomerexplorer != null)
            s = this.monomerexplorer.selected[a];
        if (!scil.Utils.isNullOrEmpty(s))
            return s;

        var set = org.helm.webeditor.Monomers.getMonomerSet(a);
        var m = set == null || this.jsd._keypresschar == null ? null : set[this.jsd._keypresschar];
        if (m != null)
            return m.id;

        if (c != null)
            return c;

        return org.helm.webeditor.Monomers.getDefaultMonomer(a);
    },

    setNodeType: function (a, biotype, elem) {
        var mon = org.helm.webeditor.Monomers.getMonomer(biotype, elem);
        if (mon != null) {
            var id = a.bio == null ? null : a.bio.id;
            a.bio = { type: biotype, id: id };
            a.elem = mon.id;
            return true;
        }

        var f = org.helm.webeditor.isAmbiguous(elem, biotype);
        if (!f) {
            elem = "(" + elem + ")";
            f = org.helm.webeditor.isAmbiguous(elem, biotype);
        }

        if (f) {
            a.elem = "?";
            a.bio = { type: biotype, id: id, ambiguity: elem };
            return true;
        }

        return false;
    },

    setNodeTypeFromGui: function (a, s) {
        var m = scil.helm.Monomers.getMonomer(a, s);
        if (m != null) {
            a.elem = s;
            return true;
        }

        var f = org.helm.webeditor.isAmbiguous(s, a.biotype());
        if (!f) {
            s = "(" + s + ")";
            f = org.helm.webeditor.isAmbiguous(s, a.biotype());
        }

        if (f) {
            a.elem = "?";
            a.bio.ambiguity = s;
            return true;
        }

        scil.Utils.alert("Unknown monomer type: " + s);
        return false;
    },

    cancelDnD: function () {
        if (this.monomerexplorer != null)
            this.monomerexplorer.dnd.cancel();
    },

    /**
    * Replace monomers
    * @function replaceMonomer
    * @param {enum} monomertype - org.helm.webeditor.HELM.BASE/SUGAR/LINKER/AA/CHEM
    * @param {string} find - the monomer name to be found
    * @param {string} replacedwith - the monomer name to be replaced with
    * @param {bool} selectedonly - indicate if seaching the selected part only
    * @returns the count replaced
    */
    replaceMonomer: function (monomertype, find, replacedwith, selectedonly) {
        var n = 0;
        var atoms = this.jsd.m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            if ((selectedonly && a.selected || !selectedonly) &&
                find == a.elem &&
                (monomertype == "" || monomertype == a.biotype())) {
                if (this.setNodeType(a, a.biotype(), replacedwith))
                    ++n;
            }
        }
        return n;
    },

    makeComplementaryStrand: function (a) {
        var chain = org.helm.webeditor.Chain.getChain(this.jsd.m, a);
        if (chain == null)
            return null;

        return chain.makeComplementaryStrand(this.jsd.m, this.jsd.bondlength);
    },

    setHelmBlobType: function (a, type) {
        if (a.biotype() == org.helm.webeditor.HELM.BLOB && a.bio.blobtype != type) {
            this.jsd.pushundo();
            a.bio.blobtype = type;
            this.jsd.refresh(true);
        }
    },

    createGroup: function (a, collapse, resetid) {
        var list = [];
        var graphics = this.jsd.m.graphics;
        for (var i = 0; i < graphics.length; ++i) {
            if (JSDraw2.Group.cast(graphics[i]) != null && graphics[i].selected)
                list.push(graphics[i]);
        }

        if (list.length > 1)
            return this.createGroup2(list, collapse);

        list = [];
        var atoms = this.jsd.m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            if (atoms[i].selected)
                list.push(atoms[i]);
        }

        if (list.length == 0)
            return false;

        var chain = null;
        var chains = org.helm.webeditor.Chain.getChains(this.jsd.m);
        if (chains.length > 0) {
            var a0 = list[0];
            for (var i = 0; i < chains.length; ++i) {
                if (scil.Utils.indexOf(chains[i].atoms, a0) >= 0) {
                    chain = chains[i];
                    break;
                }
            }
        }

        if (chain == null)
            return false;

        for (var i = 0; i < list.length; ++i) {
            if (scil.Utils.indexOf(chain.atoms, list[i]) < 0 && scil.Utils.indexOf(chain.bases, list[i]) < 0) {
                scil.Utils.alert("All monomer should be in one chain");
                return false;
            }
        }

        var openbonds = [];
        for (var i = 0; i < chain.bonds.length; ++i) {
            var b = chain.bonds[i];
            var i1 = scil.Utils.indexOf(list, b.a1);
            var i2 = scil.Utils.indexOf(list, b.a2);
            if (i1 < 0 && i2 >= 0 || i1 >= 0 && i2 < 0)
                openbonds.push({ bond: b, ai: i1 >= 0 ? 1 : 2 });
        }

        var g = this.createGroup2(list, collapse);
        if (g != null && openbonds.length > 0) {
            var sa = this.collapseGroup(g);
            org.helm.webeditor.Layout.resetIDs(sa.superatom, true);
            this.groupExpand(sa);

            for (var i = 0; i < openbonds.length; ++i) {
                var b = openbonds[i];
                if (b.ai == 1) {
                    var a1 = b.bond.a1;
                    b.bond.a1 = sa;
                    b.bond.r1 = a1._aaid + ":R" + b.bond.r1;
                }
                else {
                    var a2 = b.bond.a2;
                    b.bond.a2 = sa;
                    b.bond.r2 = a2._aaid + ":R" + b.bond.r2;
                }
            }
        }

        if (g != null && resetid)
            org.helm.webeditor.Layout.resetIDs(this.jsd.m);
        return g != null;
    },

    createGroup2: function (list, collapse) {
        var g = new JSDraw2.Group("", "helmgroup");
        g.gap = 10;
        this.jsd.m.addGraphics(g);

        for (var i = 0; i < list.length; ++i)
            list[i].group = g;

        //this.jsd.refresh(true);
        if (collapse)
            this.collapseGroup(g);
        return g;
    },

    collapseGroup: function (g, clean) {
        if (JSDraw2.Group.cast(g) == null)
            return null;

        this.jsd.m.clearFlag();
        var ret = this._collapseGroup(g);
        if (ret != null && clean)
            this.clean();

        return ret;
    },

    _collapseGroup: function (g) {
        var sa = g.a != null ? g.a : org.helm.webeditor.Interface.createAtom(this.jsd.m, new JSDraw2.Point());
        //sa.tag = "Group";
        var mol = new JSDraw2.Mol();
        sa.superatom = mol;
        sa.hidden = null;
        sa.group = g.group;
        sa.ratio = g.ratio;
        sa.tag = g.tag;
        this.setNodeType(sa, org.helm.webeditor.HELM.BLOB, "group");

        var graphics = scil.clone(this.jsd.m.graphics);
        for (var i = 0; i < graphics.length; ++i) {
            var g2 = graphics[i];
            if (g2.group == g)
                this._collapseGroup(g2);
        }

        var atoms = this.jsd.m.atoms;
        for (var i = atoms.length - 1; i >= 0; --i) {
            var a = atoms[i];
            if (a.group == g) {
                mol.atoms.push(a);
                a.f = true;
                a.selected = false;
                atoms.splice(i, 1);
            }
        }

        var apo = 0;
        var connections = [];
        var bonds = this.jsd.m.bonds;
        for (var i = bonds.length - 1; i >= 0; --i) {
            var b = bonds[i];
            if (b.a1.f != b.a2.f) {
                ++apo;
                connections.push(b);
                if (b.a1.f) {
                    b.a1.attachpoints.push(apo);
                    b.a1 = sa;
                    b.apo1 = apo;
                }
                else {
                    b.a2.attachpoints.push(apo);
                    b.a2 = sa;
                    b.apo2 = apo;
                }
            }
            else if (b.a1.f && b.a2.f) {
                mol.bonds.push(b);
                bonds.splice(i, 1);
            }
        }

        var r = org.helm.webeditor.Layout.getRect(mol.atoms);
        sa.p = r.center();

        this.jsd.m.delGraphics(g);
        //this.clean();

        //this.jsd.refresh(true);
        return sa;
    },

    groupExpand: function (a) {
        if (a == null || a.superatom == null)
            return false;

        var m = a.superatom;
        var bonds = this.jsd.m.getAllBonds(a);
        for (var i = 0; i < bonds.length; ++i) {
            var b = bonds[i];
            if (b.a1 == a) {
                var oa = this.findAtomByAP(m, b.apo1);
                if (oa != null) {
                    b.a1 = oa;
                    scil.Utils.delFromArray(oa.attachpoints, b.apo1);
                    b.apo1 = null;
                }
            }
            else {
                var oa = this.findAtomByAP(m, b.apo2);
                if (oa != null) {
                    b.a2 = oa;
                    scil.Utils.delFromArray(oa.attachpoints, b.apo2);
                    b.apo2 = null;
                }
            }
        }

        var c = m.rect().center();
        m.offset(a.p.x - c.x, a.p.y - c.y);

        this.jsd.m.mergeMol(m);

        var g = new JSDraw2.Group("", "helmgroup");
        g.gap = 10;
        this.jsd.m.addGraphics(g);
        for (var i = 0; i < m.atoms.length; ++i)
            m.atoms[i].group = g;

        g.a = a;
        g.group = a.group;
        g.ratio = a.ratio;
        g.tag = a.tag;
        a.superatom = null;
        a.hidden = true;
        //this.jsd.m.delAtom(a);
        //this.clean();

        //this.jsd.refresh(true);
        return true;
    },

    findAtomByAP: function (m, apo) {
        for (var i = 0; i < m.atoms.length; ++i) {
            var a = m.atoms[i];
            if (scil.Utils.indexOf(a.attachpoints, apo) >= 0)
                return a;
        }
        return null;
    },

    /**
    * Apply a rule
    * @function applyRule
    * @param {function} rulefun - a rule function to be called: function(plugin) {}
    */
    applyRule: function (rulefun) {
        org.helm.webeditor.RuleSet.applyRule(this, rulefun);
    },

    applyRules: function (funs) {
        org.helm.webeditor.RuleSet.applyRules(this, funs);
    },

    addNode: function (p, biotype, elem) {
        elem = org.helm.webeditor.IO.trimBracket(elem);

        var m = org.helm.webeditor.Monomers.getMonomer(biotype, elem);
        if (m == null)
            m = org.helm.webeditor.Monomers.addSmilesMonomer(biotype, elem);

        var ambiguity = null;
        if (m == null && org.helm.webeditor.isAmbiguous(elem, biotype)) {
            m = org.helm.webeditor.Monomers.getMonomer(biotype, "?");
            ambiguity = elem;
        }

        if (m == null) {
            console.warn("Unknown " + biotype + " monomer name: " + elem);
            return null;
        }

        var a = org.helm.webeditor.Interface.createAtom(this.jsd.m, p);
        this.setNodeType(a, biotype, m.id == null ? elem : m.id);
        a.bio.ambiguity = ambiguity;
        return a;
    },

    addBond: function (a1, a2, r1, r2) {
        if (a1 == null || a2 == null || a1 == a2 || r1 != "?" && !this.hasSpareR(a1, r1) || r2 != "?" && !this.hasSpareR(a2, r2))
            return null;
        //if (a1.biotype() == org.helm.webeditor.HELM.SUGAR && a2.biotype() == org.helm.webeditor.HELM.SUGAR || a1.biotype() == org.helm.webeditor.HELM.AA && a2.biotype() == org.helm.webeditor.HELM.AA) {
        //    if ((r1 == 1 || r1 == 2) && r1 == r2)
        //        return null;
        //}
        var b = org.helm.webeditor.Interface.createBond(this.jsd.m, a1, a2);
        b.r1 = r1;
        b.r2 = r2;
        return b;
    },

    addHydrogenBond: function (a1, a2) {
        if (a1 == null || a2 == null || a1 == a2)
            return null;
        var b = org.helm.webeditor.Interface.createBond(this.jsd.m, a1, a2);
        b.type = JSDraw2.BONDTYPES.UNKNOWN;
        return b;
    },

    connnectGroup: function (p1, object) {
        var object1 = this.jsd.toggle(p1);
        var object2 = object;
        var t1 = JSDraw2.Atom.cast(object1);
        var t2 = JSDraw2.Atom.cast(object2);

        var showmsg = false;
        if (t1 == null) {
            var g1 = JSDraw2.Group.cast(object1);
            if (g1 != null && !scil.Utils.isNullOrEmpty(g1.ratio)) {
                showmsg = true;
                g1.ratio = null;
            }
            t1 = this.collapseGroup(g1);
            this.groupExpand(t1);
        }

        if (t2 == null) {
            var g2 = JSDraw2.Group.cast(object2);
            if (g2 != null && !scil.Utils.isNullOrEmpty(g2.ratio)) {
                showmsg = true;
                g2.ratio = null;
            }
            t2 = this.collapseGroup(g2);
            this.groupExpand(t2);
        }

        if (showmsg)
            scil.Utils.alert("Ratios on groups are removed");

        if (t1 != null && t2 != null) {
            this.connectFragment(t1, t2);
            this.jsd.moveCenter();
            return true;
        }

        return false;
    },

    connectFragment: function (a1, a2, extendchain) {
        var b = this.jsd.m.findBond(a1, a2);
        if (b != null)
            return;

        var a = null;
        var frag = null;

        var left = a1.p.x < a2.p.x ? a1 : a2;
        if (a1.p.x > a2.p.x) {
            var t = a1;
            a1 = a2;
            a2 = t;
        }

        var delta = org.helm.webeditor.bondscale * this.jsd.bondlength;

        var bt1 = a1.biotype();
        var bt2 = a2.biotype();
        if (bt1 == org.helm.webeditor.HELM.LINKER && bt2 == org.helm.webeditor.HELM.SUGAR || bt1 == org.helm.webeditor.HELM.SUGAR && bt2 == org.helm.webeditor.HELM.LINKER || bt1 == org.helm.webeditor.HELM.SUGAR && bt2 == org.helm.webeditor.HELM.SUGAR ||
            bt1 == org.helm.webeditor.HELM.AA && bt2 == org.helm.webeditor.HELM.AA) {
            var f = false;
            if (this.hasSpareR(a1, 2) && this.hasSpareR(a2, 1)) {
                f = true;
            }
            else if (this.hasSpareR(a2, 2) && this.hasSpareR(a1, 1)) {
                var t = a1;
                a1 = a2;
                a2 = t;

                f = true;
            }

            if (f) {
                frag = this.jsd.getFragment(a2);
                if (bt1 == org.helm.webeditor.HELM.AA) {
                    b = this.addBond(a1, a2, 2, 1);
                }
                else {
                    if (bt1 != bt2 || !this.needLinker()) {
                        b = this.addBond(a1, a2, 2, 1);
                    }
                    else {
                        a = this.addNode(org.helm.webeditor.Interface.createPoint(left.p.x + delta, left.p.y), org.helm.webeditor.HELM.LINKER, this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER));
                        b = this.addBond(a1, a, 2, 1);
                        if (b != null)
                            b = this.addBond(a, a2, 2, 1);
                    }
                }

                this.finishConnect(extendchain, b, a, a1, a2, frag, delta);
                return;
            }
        }
        else if (bt1 == org.helm.webeditor.HELM.SUGAR && bt2 == org.helm.webeditor.HELM.BASE || bt2 == org.helm.webeditor.HELM.SUGAR && bt1 == org.helm.webeditor.HELM.BASE) {
            if (bt2 == org.helm.webeditor.HELM.SUGAR) {
                var t = a1;
                a1 = a2;
                a2 = t;
            }
            var b = this.addBond(a1, a2, 3, 1);
            if (b != null) {
                a2.p = org.helm.webeditor.Interface.createPoint(a1.p.x, a1.p.y + Math.abs(delta));
                this.finishConnect(false, b, null, b.a1);
            }
            return;
        }

        var rs1 = this.getSpareRs(a1);
        var rs2 = this.getSpareRs(a2);
        if ((rs1 == null || rs2 == null)) {
            if (this.canPair(a1, a2) && this.jsd.m.findBond(a1, a2) == null) {
                // hydrogen bond
                org.helm.webeditor.Interface.createBond(this.jsd.m, a1, a2, JSDraw2.BONDTYPES.UNKNOWN);
                this.finishConnect(extendchain);
            }
            else {
                var s = "";
                if (rs1 == null && rs2 == null)
                    s = "Both monomers don't";
                else if (rs1 == null)
                    s = "Monomer, " + a1.elem + (a1.bio.id == null ? "" : a1.bio.id) + ", doesn't";
                else if (rs2 == null)
                    s = "Monomer, " + a2.elem + (a2.bio.id == null ? "" : a2.bio.id) + ", doesn't";
                scil.Utils.alert(s + " have any connecting point available");
                this.finishConnect(extendchain);
            }
            return;
        }

        if (rs1.length <= 1 && rs2.length <= 1) {
            if (bt1 == org.helm.webeditor.HELM.LINKER)
                bt1 = org.helm.webeditor.HELM.SUGAR;
            if (bt2 == org.helm.webeditor.HELM.LINKER)
                bt2 = org.helm.webeditor.HELM.SUGAR;

            // https://github.com/PistoiaHELM/HELMWebEditor/issues/101
            // prevent head-to-head and tail-to-tail connection
            if (!org.helm.webeditor.allowHeadToHeadConnection) {
                if (bt1 == bt2 && (bt1 == org.helm.webeditor.HELM.SUGAR || bt1 == org.helm.webeditor.HELM.AA) && rs1[0] == rs2[0] && (rs1[0] == 1 || rs1[0] == 2)) {
                    scil.Utils.alert("head-to-head / tail-to-tail connection is not allowed");
                    return;
                }
            }

            frag = this.jsd.getFragment(a2);
            b = this.addBond(a1, a2, rs1[0], rs2[0]);
        }
        else {
            if (extendchain)
                this.jsd.refresh();

            var me = this;
            this.chooseRs(rs1, rs2, a1, a2, function (r1, r2) {
                frag = me.jsd.getFragment(a2);
                b = me.addBond(a1, a2, r1, r2);
                me.finishConnect(extendchain, b, a1, a1, a2, frag, delta);
            });
            return;
        }

        this.finishConnect(extendchain, b, a, a1, a2, frag, delta);
    },

    canPair: function (a1, a2) {
        if (a1.biotype() == org.helm.webeditor.HELM.BASE && a2.biotype() == org.helm.webeditor.HELM.BASE) {
            var c1 = a1.elem;
            var c2 = a2.elem;
            return c1 == "A" && (c2 == "T" || c2 == "U") || (c1 == "T" || c1 == "U") && c2 == "A" ||
                c1 == "G" && c2 == "C" || c1 == "C" && c2 == "G";
        }
        return false;
    },

    needLinker: function () {
        var linker = this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER);
        return linker != "null";
    },

    finishConnect: function (extendchain, b, a, a1, a2, frag, delta) {
        this.clean();
        this.jsd.refresh(extendchain || b != null);
    },

    chooseRs: function (rs1, rs2, a1, a2, callback) {
        if (this.chooseRDlg == null) {
            var me = this;
            var fields = {
                s1: { label: "Monomer 1", type: "jsdraw", width: 240, height: 150, viewonly: true, style: { border: "solid 1px gray"} },
                r1: { type: "select", width: 120 },
                g: { type: "div" },
                s2: { label: "Monomer 2", type: "jsdraw", width: 240, height: 150, viewonly: true, style: { border: "solid 1px gray"} },
                r2: { type: "select", width: 120 }
            };
            this.chooseRDlg = scil.Form.createDlgForm("Choose Connecting Points", fields, { label: "OK", width: 240, onclick: function () { me.chooseRs2(); } });
        }

        this.chooseRDlg.callback = callback;
        this.chooseRDlg.show2({ owner: this.jsd });
        this._listRs(this.chooseRDlg.form.fields.r1, rs1, 2);
        this._listRs(this.chooseRDlg.form.fields.r2, rs2, 1);

        this.chooseRDlg.form.fields.r1.disabled = rs1.length <= 1;
        this.chooseRDlg.form.fields.r2.disabled = rs2.length <= 1;

        var m1 = org.helm.webeditor.Monomers.getMonomer(a1);
        var m2 = org.helm.webeditor.Monomers.getMonomer(a2);
        this.chooseRDlg.form.fields.s1.jsd.setMolfile(org.helm.webeditor.Monomers.getMolfile(m1));
        this.chooseRDlg.form.fields.s2.jsd.setMolfile(org.helm.webeditor.Monomers.getMolfile(m2));

        var tr1 = scil.Utils.getParent(this.chooseRDlg.form.fields.s1, "TR");
        var tr2 = scil.Utils.getParent(this.chooseRDlg.form.fields.s2, "TR");
        tr1.childNodes[0].innerHTML = a1.elem + (a1.bio == null || a1.bio.id == null ? "" : a1.bio.id);
        tr2.childNodes[0].innerHTML = a2.elem + (a2.bio == null || a2.bio.id == null ? "" : a2.bio.id);

        this.chooseRDlg.rs1 = rs1;
        this.chooseRDlg.rs2 = rs2;
    },

    chooseRs2: function () {
        var d = this.chooseRDlg.form.getData();
        if (scil.Utils.isNullOrEmpty(d.r1) && this.chooseRDlg.rs1.length > 0 || scil.Utils.isNullOrEmpty(d.r2) && this.chooseRDlg.rs2.length > 0) {
            scil.Utils.alert("Please select Rs for both Nodes");
            return;
        }

        this.chooseRDlg.hide();
        this.chooseRDlg.callback(d.r1 == null || d.r1 == "?" ? d.r1 : parseInt(d.r1), d.r2 == null || d.r2 == "?" ? d.r2 : parseInt(d.r2));
    },

    _listRs: function (sel, list, v) {
        var ss = {};
        for (var i = 0; i < list.length; ++i)
            ss[list[i] + ""] = list[i] == "?" ? "?" : ("R" + list[i]);
        scil.Utils.listOptions(sel, ss, v == null ? null : (v + ""), true, false);
    },

    changeMonomer: function (a, cloned) {
        var s = this.getDefaultNodeType(a.biotype());
        if (!scil.Utils.isNullOrEmpty(s) && a.elem != s && s != "null" && a.biotype() != org.helm.webeditor.HELM.BLOB) {
            this.jsd.pushundo(cloned);
            this.setNodeType(a, a.biotype(), s);
            this.jsd.refresh(true);
        }
        else {
            scil.Utils.beep();
        }
    },

    extendChain: function (a1, cmd, p1, p2, cloned) {
        var rs = [];
        var rgroups = this.getSpareRs(a1, rs);
        if (rgroups == null) {
            scil.Utils.alert("No connecting points available");
            this.jsd.redraw();
            return;
        }

        var delta = p2.x > p1.x ? org.helm.webeditor.bondscale * this.jsd.bondlength : -org.helm.webeditor.bondscale * this.jsd.bondlength;
        var p = org.helm.webeditor.Interface.createPoint(a1.p.x + delta, a1.p.y);

        var a2 = null;
        var r1 = null;
        var r2 = null;
        var bond = null;
        if (cmd == "helm_chem") {
            if (Math.abs(p2.y - p1.y) / Math.abs(p2.x - p1.x) > 5)
                p = org.helm.webeditor.Interface.createPoint(a1.p.x, a1.p.y + Math.abs(delta) * (p2.y > p1.y ? 1 : -1));
            a2 = this.addNode(p, org.helm.webeditor.HELM.CHEM, this.getDefaultNodeType(org.helm.webeditor.HELM.CHEM));
            if (a2 != null) {
                this.connectFragment(a1, a2, true);
                return;
            }
        }
        else if (cmd == "helm_aa") {
            if (Math.abs(p2.y - p1.y) / Math.abs(p2.x - p1.x) > 5)
                p = org.helm.webeditor.Interface.createPoint(a1.p.x, a1.p.y + Math.abs(delta) * (p2.y > p1.y ? 1 : -1));
            a2 = this.addNode(p, org.helm.webeditor.HELM.AA, this.getDefaultNodeType(org.helm.webeditor.HELM.AA));
        }
        else if (cmd == "helm_linker") {
            a2 = this.addNode(p, org.helm.webeditor.HELM.LINKER, this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER));
        }
        else if (cmd == "helm_sugar") {
            a2 = this.addNode(p, org.helm.webeditor.HELM.SUGAR, this.getDefaultNodeType(org.helm.webeditor.HELM.SUGAR));
        }
        else if (cmd == "helm_base") {
            if (a1.biotype() == org.helm.webeditor.HELM.SUGAR && this.hasSpareR(a1, 3)) {
                r1 = 3;
                r2 = 1;
                p = org.helm.webeditor.Interface.createPoint(a1.p.x, a1.p.y + Math.abs(delta));
                a2 = this.addNode(p, org.helm.webeditor.HELM.BASE, this.getDefaultNodeType(org.helm.webeditor.HELM.BASE));
            }
        }
        else if (cmd == "helm_nucleotide" || cmd == "helm_sugar") {
            if (Math.abs(p2.y - p1.y) / Math.abs(p2.x - p1.x) > 5) {
                // drag vertically to add base
                if (a1.biotype() == org.helm.webeditor.HELM.SUGAR && rs[3] == true) {
                    r1 = 3;
                    r2 = 1;
                    p = org.helm.webeditor.Interface.createPoint(a1.p.x, a1.p.y + Math.abs(delta));
                    a2 = this.addNode(p, org.helm.webeditor.HELM.BASE, this.getDefaultNodeType(org.helm.webeditor.HELM.BASE));
                }
            }
            else {
                if (rs[1] == true || rs[2] == true) {
                    var m = this.getDefaultNodeType(org.helm.webeditor.HELM.SUGAR);
                    var e = this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER);
                    var linker = null;
                    var sugar = null;

                    if (delta < 0) {
                        if (rs[1])
                            r1 = 1;
                        else
                            r1 = 2;
                    }
                    else {
                        if (rs[2])
                            r1 = 2;
                        else
                            r1 = 1;
                    }
                    r2 = r1 == 1 ? 2 : 1;

                    if (r1 == 1) {
                        if (e != "null") {
                            linker = this.addNode(p.clone(), org.helm.webeditor.HELM.LINKER, e);
                            p.x += delta;
                        }
                        sugar = this.addNode(p.clone(), org.helm.webeditor.HELM.SUGAR, m);

                        if (linker != null) {
                            bond = this.addBond(a1, linker, r1, r2);
                            this.addBond(linker, sugar, r1, r2);
                        }
                        else {
                            bond = this.addBond(a1, sugar, r1, r2);
                        }
                    }
                    else {
                        sugar = this.addNode(p.clone(), org.helm.webeditor.HELM.SUGAR, m);
                        p.x += delta;
                        if (e != "null")
                            linker = this.addNode(p.clone(), org.helm.webeditor.HELM.LINKER, e);

                        if (linker != null) {
                            bond = this.addBond(a1, sugar, r1, r2);
                            this.addBond(sugar, linker, r1, r2);
                        }
                        else {
                            bond = this.addBond(a1, sugar, r1, r2);
                        }
                    }

                    var base = null;
                    if (cmd == "helm_nucleotide" && org.helm.webeditor.Monomers.hasR(org.helm.webeditor.HELM.SUGAR, m, "R3")) {
                        base = this.addNode(sugar.p.clone(), org.helm.webeditor.HELM.BASE, this.getDefaultNodeType(org.helm.webeditor.HELM.BASE));
                        this.addBond(sugar, base, 3, 1);

                        var leftR = bond.a1.p.x < bond.a2.p.x ? bond.r1 : bond.r2;
                        if (leftR == 1) // reversed
                            base.p.y -= Math.abs(delta);
                        else
                            base.p.y += Math.abs(delta);
                    }

                    this.jsd.pushundo(cloned);
                    this.finishConnect(false, bond, null, a1);
                }
            }
        }

        if (a2 != null) {
            if (r1 == null || r2 == null) {
                if (this.hasSpareR(a1, 2) && !this.hasSpareR(a1, 1)) {
                    r1 = 2;
                    r2 = 1;
                }
                else if (this.hasSpareR(a1, 1) && !this.hasSpareR(a1, 2)) {
                    r1 = 1;
                    r2 = 2;
                }
                else {
                    r1 = delta > 0 ? 2 : 1;
                    r2 = r1 == 2 ? 1 : 2;
                }
            }
            bond = this.addBond(a1, a2, r1, r2);
        }

        if (bond != null) {
            this.jsd.pushundo(cloned);
            this.finishConnect(false, bond, null, a1);
        }
        else {
            this.jsd.refresh();
        }
    },

    /**
    * Get HELM
    * @function getHelm
    * @param {bool} highlightselection - internal use only, using null always
    * @returns the HELM as string
    */
    getHelm: function (highlightselection) {
        return org.helm.webeditor.IO.getHelm(this.jsd.m, highlightselection);
    },

    /**
    * Get sequence of natuaral analogue
    * @function getSequence
    * @param {bool} highlightselection - internal use only, using null always
    * @returns the sequence as a string
    */
    getSequence: function (highlightselection) {
        return org.helm.webeditor.IO.getSequence(this.jsd.m, highlightselection);
    },

    /**
    * Get XHELM
    * @function getXHelm
    * @returns XHELM as a string
    */
    getXHelm: function () {
        return org.helm.webeditor.IO.getXHelm(this.jsd.m);
    },

    /**
    * Set HELM
    * @function getSequence
    * @param {string} s - The HELM string
    * @param {string} renamedmonomers - internal use only, using null always
    */
    setHelm: function (s, renamedmonomers) {
        this.jsd.clear();

        var n = 0;
        try {
            if (!scil.Utils.isNullOrEmpty(s))
                n = org.helm.webeditor.IO.read(this, s, "HELM", renamedmonomers);
        }
        catch (e) {
            this.jsd.clear();
            return;
        }

        if (n > 0) {
            this.clean();
            this.jsd.fitToWindow();
            this.jsd.refresh();
        }
    },

    /**
    * Set XHELM 
    * @function setXHelm
    * @param {string} s - the xhelm string
    */
    setXHelm: function (s) {
        var doc = typeof (s) == "string" ? scil.Utils.parseXml(s) : s;
        if (doc == null)
            return false;

        var es = doc.getElementsByTagName("HelmNotation");
        if (es == null || es.length == 0)
            return false;

        var s = scil.Utils.getInnerText(es[0]);

        var list = doc.getElementsByTagName("Monomers");
        if (list == null || list.length == 0) {
            this.setHelm(s);
            return;
        }

        var me = this;
        org.helm.webeditor.monomers.loadMonomers(list[0], function (renamed) {
            if (me.monomerexplorer != null)
                me.monomerexplorer.reloadTabs();
            me.setHelm(s, renamed);
        });
    },

    isXHelm: function (doc) {
        var ret = doc == null ? null : doc.getElementsByTagName("Xhelm");
        return ret != null && ret.length == 1;
    },

    /**
    * Show Importing Sequence dialog
    * @function showImportDlg
    */
    showImportDlg: function () {
        if (this.inputSeqDlg == null) {
            var fields = {
                type: { label: "Sequence Type", type: "select", items: ["HELM", "Peptide", "RNA"] },
                sequence: { label: "Sequence", type: "textarea", width: 800, height: 50 },
                separator: { label: "Separator", width: 100, str: "(Legacy sequences with dedicated separator, e.g. dC.dA.xE)" }
            };

            var me = this;
            this.inputSeqDlg = scil.Form.createDlgForm("Import Sequence", fields, [
                { label: "Import", onclick: function () { me.importSequence(false); } },
                { label: "Append", onclick: function () { me.importSequence(true); } }
            ]);
        }

        this.inputSeqDlg.show2({ owner: this.jsd });
    },

    importSequence: function (append) {
        var data = this.inputSeqDlg.form.getData();
        if (this.setSequence(data.sequence, data.type, null, null, append, data.separator))
            this.inputSeqDlg.hide();
    },

    /**
    * Set a sequence (natural analogue sequence, HELM,)
    * @function setSequence
    * @param {string} seq - the input sequence
    * @param {string} format - input format: HELM, RNA, Peptide, or null
    * @param {string} sugar - the sugar for RNA
    * @param {string} linker - the linker for RNA
    * @param {bool} append - set the sequence in appending mode or overwriting mode
    */
    setSequence: function (seq, format, sugar, linker, append, separator) {
        var seq = scil.Utils.trim(seq);
        if (/^[a-z]+$/.test(seq))
            seq = seq.toUpperCase();

        var n = 0;
        var cloned = this.jsd.clone();
        this.jsd.clear();
        try {
            n = org.helm.webeditor.IO.read(this, seq, format, null, sugar, linker, separator);
        }
        catch (e) {
            this.jsd.restoreClone(cloned);
            var s = e.message == null ? e : e.message;
            if (!scil.Utils.isNullOrEmpty(s))
                scil.Utils.alert("Error: " + s);
            return false;
        }

        if (n > 0) {
            this.jsd.pushundo(cloned);

            this.clean();

            if (append) {
                var m = cloned.mol.clone();
                var rect = m.rect();
                var r2 = this.jsd.m.rect();
                if (r2 != null && rect != null)
                    this.jsd.m.offset(rect.center().x - r2.center().x, rect.bottom() + this.jsd.bondlength * 4 - r2.bottom());
                m.mergeMol(this.jsd.m);
                this.jsd.m = m;
            }

            this.jsd.fitToWindow();
            this.jsd.refresh(true);
        }
        return true;
    },

    /**
    * Clean the layout
    * @function clean
    * @param {JSDraw2.Atom} a - the start monomer.  Use null to clean all
    * @param {bool} redraw - indicate if redrawing the structure after cleaning
    */
    clean: function (a, redraw) {
        if (redraw)
            this.jsd.pushundo();

        org.helm.webeditor.Layout.clean(this.jsd, a);
        if (redraw) {
            this.jsd.moveCenter();
            this.jsd.refresh(true);
        }
    },

    /**
    * Reset monomer IDs 
    * @function resetIDs
    */
    resetIDs: function () {
        org.helm.webeditor.Layout.resetIDs(this.jsd.m);
    },

    dropMonomer: function (type, id, e) {
        var p = this.jsd.eventPoint(e);
        if (p.x <= 0 || p.y <= 0 || p.x >= this.jsd.dimension.x || p.y >= this.jsd.dimension.y || id == "null")
            return false;

        var f = false;
        if (this.jsd.curObject == null) {
            // create new monomer
            var cmd = type == "nucleotide" ? "helm_nucleotide" : scil.helm.symbolCase(type);
            if (this.isHelmCmd(cmd)) {
                p.offset(this.jsd.bondlength * 0.4, this.jsd.bondlength * 0.4);
                this.jsd.pushundo();
                var a = org.helm.webeditor.Interface.createAtom(this.jsd.m, p);
                this.createIsolatedMonomer(cmd, a);
                f = true;
            }
        }
        else {
            // modify the target monomer
            var set = org.helm.webeditor.Monomers.getMonomerSet(type);
            if (org.helm.webeditor.isAmbiguous(id, type)) {
                var a = org.helm.webeditor.Interface.getCurrentAtom(this.jsd);
                if (a == null || !org.helm.webeditor.isHelmNode(a) || a.biotype() != type || a.elem == id)
                    return false;
            }
            else {
                if (set == null || set[scil.helm.symbolCase(id)] == null)
                    return false;

                id = set[scil.helm.symbolCase(id)].id;
                var a = org.helm.webeditor.Interface.getCurrentAtom(this.jsd);
                if (a == null || !org.helm.webeditor.isHelmNode(a) || a.biotype() != type || a.elem == id)
                    return false;
            }
            this.jsd.pushundo();
            this.setNodeType(a, a.biotype(), id);
            f = true;
        }

        if (f)
            this.jsd.refresh(true);
        return f;
    },

    showFindReplaceDlg: function () {
        if (this.findDlg == null) {
            var fields = {
                finding: { label: "Find", width: 400, str: "<div>(Monomer symbol or position)</div>" },
                monomertype: { label: "Monomer Type", type: "select", sort: false, items: org.helm.webeditor.monomerTypeList() },
                replacewith: { label: "Replace With", width: 400 },
                selectedonly: { label: "Scope", type: "checkbox", str: "Search Selected Only" }
            };

            var me = this;
            this.findDlg = scil.Form.createDlgForm("Find and Replace", fields, [
                { label: "Find", onclick: function () { me.showFindReplaceDlg2("find"); } },
                { label: "Find All", onclick: function () { me.showFindReplaceDlg2("findall"); } },
                { label: "Replace All", onclick: function () { me.showFindReplaceDlg2("replaceall"); } }
            ])
        }

        this.findDlg.show2({ owner: this.jsd });
    },

    showFindReplaceDlg2: function (action) {
        var data = this.findDlg.form.getData();
        if (scil.Utils.isNullOrEmpty(data.finding) || action == "replaceall" && scil.Utils.isNullOrEmpty(data.replacewith)) {
            scil.Utils.alert("Find and Replace With cannot be blank");
            return;
        }

        if (action == "find")
            this.find(data.finding, false, data.monomertype, data.selectedonly);
        else if (action == "findall")
            this.find(data.finding, true, data.monomertype, data.selectedonly);
        else if (action == "replaceall")
            this.replaceAll(data.finding, data.replacewith, data.monomertype, data.selectedonly);
    },

    getSelectedAtoms: function () {
        var ret = [];
        var atoms = this.jsd.m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            if (atoms[i].selected)
                ret.push(atoms[i]);
        }
        return ret;
    },

    find: function (a, findall, monomertype, selectedonly) {
        var atoms = selectedonly ? this.getSelectedAtoms() : this.jsd.m.atoms;
        this.jsd.m.setSelected(false);

        var n = 0;
        var atom = null;
        if (/^[0-9]+$/.test(a)) {
            var aaid = parseInt(a);
            for (var i = 0; i < atoms.length; ++i) {
                if (atoms[i].bio != null && aaid == atoms[i].bio.id && (scil.Utils.isNullOrEmpty(monomertype) || monomertype == atoms[i].biotype())) {
                    ++n;
                    atoms[i].selected = true;
                    atom = atoms[i];
                    break;
                }
            }
        }
        else {
            for (var i = 0; i < atoms.length; ++i) {
                if (a == atoms[i].elem && (monomertype == "" || monomertype == atoms[i].biotype())) {
                    ++n;
                    atoms[i].selected = true;
                    if (!findall) {
                        atom = atoms[i];
                        break;
                    }
                }
            }
        }

        if (findall) {
            scil.Utils.alert(n + " node(s) found");
        }
        else {
            if (n == 0) {
                scil.Utils.alert("Cannot find " + a);
            }
            else {
                org.helm.webeditor.Interface.scaleCanvas(this.jsd);
                var dx = this.jsd.dimension.x / 2 - atom.p.x;
                var dy = this.jsd.dimension.y / 2 - atom.p.y;
                this.jsd.m.offset(dx, dy);
            }
        }

        if (n > 0)
            this.jsd.redraw();
    },

    replaceAll: function (a, a2, monomertype, selectedonly) {
        var n = 0;
        var cloned = this.jsd.clone();
        if (/^[0-9]+$/.test(a)) {
            var aaid = parseInt(a);
            var atoms = selectedonly ? this.getSelectedAtoms() : this.jsd.m.atoms;
            for (var i = 0; i < atoms.length; ++i) {
                if (atoms[i].bio != null && aaid == atoms[i].bio.id && (monomertype == "" || monomertype == atoms[i].biotype())) {
                    if (this.setNodeType(atoms[i], atoms[i].biotype(), a2))
                        ++n;
                    break;
                }
            }
        }
        else {
            n = this.replaceMonomer(monomertype, a, a2, selectedonly);
        }

        if (n > 0) {
            this.jsd.pushundo(cloned);
            this.jsd.refresh(true);
        }

        scil.Utils.alert(n + " node(s) replaced");
    },

    dblclickMomonor: function (type, monomer) {
        if (monomer == "null")
            return;

        var list = [];
        var atoms = this.jsd.m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            if (atoms[i].selected && atoms[i].biotype() == type && atoms[i].elem != monomer)
                list.push(atoms[i]);
        }

        if (list.length > 0) {
            this.jsd.pushundo();
            for (var i = 0; i < list.length; ++i)
                this.setNodeType(list[i], list[i].biotype(), monomer);
            this.jsd.refresh(true);
        }

        return list.length;
    },

    isHelmCmd: function (cmd) {
        return cmd == "helm_nucleotide" || cmd == "helm_base" || cmd == "helm_sugar" || cmd == "helm_chem" ||
            cmd == "helm_aa" || cmd == "helm_linker" || cmd == "helm_blob";
    },

    createIsolatedMonomer: function (cmd, a) {
        if (cmd == "helm_nucleotide") {
            var s = this.monomerexplorer == null ? null : scil.Utils.getInnerText(this.monomerexplorer.lastdiv);
            if (s == "*") {
                this.setNodeType(a, org.helm.webeditor.HELM.NUCLEOTIDE, "*");
                return true;
            }

            var m = this.getDefaultNodeType(org.helm.webeditor.HELM.SUGAR);
            this.setNodeType(a, org.helm.webeditor.HELM.SUGAR, m);

            if (org.helm.webeditor.Monomers.hasR(org.helm.webeditor.HELM.SUGAR, m, "R3")) {
                var a3 = this.addNode(org.helm.webeditor.Interface.createPoint(a.p.x, a.p.y + this.jsd.bondlength * org.helm.webeditor.bondscale), org.helm.webeditor.HELM.BASE, this.getDefaultNodeType(org.helm.webeditor.HELM.BASE));
                this.addBond(a, a3, 3, 1);
            }

            var linker = this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER);
            if (linker == null || linker == "null")
                return;

            var a2 = this.addNode(org.helm.webeditor.Interface.createPoint(a.p.x + this.jsd.bondlength * org.helm.webeditor.bondscale, a.p.y), org.helm.webeditor.HELM.LINKER, linker);
            this.addBond(a, a2, 2, 1);
        }
        else if (cmd == "helm_base") {
            this.setNodeType(a, org.helm.webeditor.HELM.BASE, this.getDefaultNodeType(org.helm.webeditor.HELM.BASE));
        }
        else if (cmd == "helm_sugar") {
            this.setNodeType(a, org.helm.webeditor.HELM.SUGAR, this.getDefaultNodeType(org.helm.webeditor.HELM.SUGAR));
        }
        else if (cmd == "helm_linker") {
            this.setNodeType(a, org.helm.webeditor.HELM.LINKER, this.getDefaultNodeType(org.helm.webeditor.HELM.LINKER));
        }
        else if (cmd == "helm_aa") {
            this.setNodeType(a, org.helm.webeditor.HELM.AA, this.getDefaultNodeType(org.helm.webeditor.HELM.AA));
        }
        else if (cmd == "helm_chem") {
            this.setNodeType(a, org.helm.webeditor.HELM.CHEM, this.getDefaultNodeType(org.helm.webeditor.HELM.CHEM));
        }
        else if (cmd == "helm_blob") {
            this.setNodeType(a, org.helm.webeditor.HELM.BLOB, this.getDefaultNodeType(org.helm.webeditor.HELM.BLOB));
            this.setHelmBlobType(a, org.helm.webeditor.blobtypes[0]);
        }
        else {
            return false;
        }

        return true;
    }
});
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* Chain class
* @class org.helm.webeditor.Chain
*/
org.helm.webeditor.Chain = scil.extend(scil._base, {
    /**
    * Contructor of Chain
    * @constructor Chain
    */
    constructor: function (id) {
        this.id = id;
        this.bonds = [];
        this.atoms = [];
        this.basebonds = [];
        this.bases = [];
    },

    /**
    * Get the complementary of a given node
    * @function getComplementary
    * @param {Atom} a - the input Node
    */
    getComplementary: function (a) {
        var m = org.helm.webeditor.Monomers.getMonomer(a);
        switch (m.na) {
            case "A":
                return "T";
            case "G":
                return "C";
            case "T":
            case "U":
                return "A";
            case "C":
                return "G";
            default:
                return "U";
        }
    },

    /**
    * Create the complementary strand (internal use)
    * @function makeComplementaryStrand
    */
    makeComplementaryStrand: function (m, bondlength) {
        var ret = new org.helm.webeditor.Chain(null);

        var lasta2 = null;
        var lastsugar = null;
        var d = bondlength * org.helm.webeditor.bondscale;
        for (var i = 0; i < this.atoms.length; ++i) {
            var a = this.atoms[i];
            var b = this.bases[i];

            var a2 = a.clone();
            ret.atoms.splice(0, 0, a2);
            a2.p.y += 3 * d;
            a2.bio.annotation = null;
            m.addAtom(a2);
            if (b != null) {
                var b2 = b.clone();
                ret.bases.splice(0, 0, b2);
                b2.p.y += d;
                b2.elem = this.getComplementary(b);
                m.addAtom(b2);

                var bond = new JSDraw2.Bond(a2, b2);
                ret.basebonds.splice(0, 0, bond);
                bond.r1 = 3;
                bond.r2 = 1;
                m.addBond(bond);

                bond = new JSDraw2.Bond(b, b2, JSDraw2.BONDTYPES.UNKNOWN);
                m.addBond(bond);
            }
            else {
                ret.bases.splice(0, 0, null);
            }

            if (lasta2 != null) {
                var bond = new JSDraw2.Bond(lasta2, a2);
                ret.bonds.splice(0, 0, bond);
                bond.r1 = 1;
                bond.r2 = 2;
                m.addBond(bond);
            }

            lasta2 = a2;
            if (a2.biotype() == org.helm.webeditor.HELM.SUGAR) {
                lastsugar = a2;
                a2.elem = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.SUGAR);
            }
            else if (a2.biotype() == org.helm.webeditor.HELM.LINKER) {
                a2.elem = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.LINKER);
            }
        }

        if (lastsugar != null)
            lastsugar.bio.annotation = "5'";

        ret.resetIDs();
        return ret;
    },

    /**
    * Get separated polymer from this fragment(chain) (internal use)
    * @function _getPolymers
    */
    _getPolymers: function () {
        var ret = [];

        var polymer = null;
        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            var biotype = a.biotype();
            if (biotype == org.helm.webeditor.HELM.AA) {
                if (polymer != null && polymer.type != "Peptide")
                    polymer = null;

                if (polymer == null) {
                    polymer = { type: "Peptide", atoms: [] };
                    ret.push(polymer);
                }
                polymer.atoms.push(a);
            }
            else if (biotype == org.helm.webeditor.HELM.SUGAR || biotype == org.helm.webeditor.HELM.HELM_LINKER) {
                if (polymer != null && polymer.type != "RNA")
                    polymer = null;

                if (biotype == org.helm.webeditor.HELM.SUGAR) {
                    var b = this.bases[i];
                    if (b != null) {
                        if (polymer == null) {
                            polymer = { type: "RNA", atoms: [] };
                            ret.push(polymer);
                        }
                        polymer.atoms.push(b);
                    }
                }
            }
            else {
                polymer = null;
            }
        }

        return ret;
    },

    /**
    * Get whole chain as an expanded molecule (internal use)
    * @function getMol
    */
    getMol: function (a, plugin) {
        var mon = org.helm.webeditor.Monomers.getMonomer(a);
        var molfile = org.helm.webeditor.monomers.getMolfile(mon);
        var m = org.helm.webeditor.Interface.createMol(molfile);

        if (plugin != null) {
            for (var r in mon.at) {
                if (plugin.hasSpareR(a, r))
                    org.helm.webeditor.MolViewer.capRGroup(m, r, mon);
            }
        }

        for (var i = 0; i < m.atoms.length; ++i)
            m.atoms[i]._helmgroup = a;

        return m;
    },

    /**
    * Expaned backbond (internal use)
    * @function _expandBackbone
    */
    _expandBackbone: function (mol, plugin) {
        var m1 = null;
        var m2 = null;
        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            var b = this.bases[i];

            a.f = true;
            m2 = this.getMol(a, plugin);

            if (b != null) {
                var m3 = this.getMol(b, plugin);
                org.helm.webeditor.MolViewer.mergeMol(m2, "R3", m3, "R1", a, b);
            }

            if (m1 == null) {
                m1 = m2;
            }
            else {
                var r1, r2;
                var a1, a2;
                var bond = this.bonds[i - 1];
                if (bond.a2 == a) {
                    r2 = bond.r2;
                    r1 = bond.r1;
                    a1 = bond.a1;
                    a2 = bond.a2;
                }
                else {
                    r2 = bond.r1;
                    r1 = bond.r2;
                    a1 = bond.a2;
                    a2 = bond.a1;
                }

                org.helm.webeditor.MolViewer.mergeMol(m1, "R" + r1, m2, "R" + r2, a1, a2);
            }
        }

        if (this.isCircle()) {
            var bond = this.bonds[n - 1];

            // circle
            var t = org.helm.webeditor.MolViewer.findR(m1, "R" + bond.r1, bond.a1);
            var s = org.helm.webeditor.MolViewer.findR(m1, "R" + bond.r2, bond.a2);
            if (t != null && s != null) {
                m1.atoms.splice(scil.Utils.indexOf(m1.atoms, t.a1), 1);
                m1.atoms.splice(scil.Utils.indexOf(m1.atoms, s.a1), 1);
                m1.bonds.splice(scil.Utils.indexOf(m1.bonds, s.b), 1);

                if (t.b.a1 == t.a1)
                    t.b.a1 = s.a0;
                else
                    t.b.a2 = s.a0;
            }
        }

        mol.atoms = mol.atoms.concat(m1.atoms);
        mol.bonds = mol.bonds.concat(m1.bonds);
    },

    /**
    * Connect branches to the backbone (internal use)
    * @function _connectBranches
    */
    _connectBranches: function (m, plugin, branches) {
        if (branches == null || branches.bonds == null)
            return;

        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            this._connectBranches2(m, a, branches, plugin);

            var b = this.bases[i];
            if (b != null)
                this._connectBranches2(m, b, branches, plugin);
        }
    },

    /**
    * Inner loop of connecting branches (internal use)
    * @function _connectBranches2
    */
    _connectBranches2: function (m, a, branches, plugin) {
        var r1 = null;
        var r2 = null;
        var a1 = null;
        var a2 = null;
        for (var i = 0; i < branches.bonds.length; ++i) {
            var b = branches.bonds[i];
            if (b == null)
                continue;

            if (b.a1 == a && !b.a2.f) {
                r1 = b.r1;
                r2 = b.r2;
                a1 = b.a1;
                a2 = b.a2;
            }
            else if (b.a2 == a && !b.a1.f) {
                r1 = b.r2;
                r2 = b.r1;
                a1 = b.a2;
                a2 = b.a1;
            }

            if (a2 != null) {
                b.f = true;
                break;
            }
        }

        if (a2 == null)
            return;

        var m2 = this.getMol(a2, plugin);
        org.helm.webeditor.MolViewer.mergeMol(m, "R" + r1, m2, "R" + r2, a, a2);
    },

    /**
    * Check if this chain is a circle
    * @function isCircle
    */
    isCircle: function () {
        return this.atoms.length >= 3 && this.atoms[0] == this.atoms[this.atoms.length - 1];
    },

    /**
    * Reset chain monomer IDs
    * @function resetIDs
    */
    resetIDs: function (resetaaid) {
        var aaid = 0;
        var _aaid = 0;
        var baseid = 0;

        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            var biotype = a.biotype();
            if (biotype == org.helm.webeditor.HELM.AA) {
                a._aaid = ++_aaid;
                a.bio.id = ++aaid;
                if (aaid == 1)
                    a.bio.annotation = "n";
                else
                    a.bio.annotation = null;
                baseid = 0;
            }
            else if (biotype == org.helm.webeditor.HELM.SUGAR || biotype == org.helm.webeditor.HELM.LINKER) {
                if (biotype == org.helm.webeditor.HELM.SUGAR && this.bases[i] != null) {
                    a._aaid = ++_aaid;
                    this.bases[i]._aaid = ++_aaid;

                    this.bases[i].bio.id = ++baseid;
                    if (baseid == 1) {
                        if (a.bio.annotation != "5'ss" && a.bio.annotation != "5'as")
                            a.bio.annotation = "5'";
                    }
                    else {
                        a.bio.annotation = null;
                    }
                }
                aaid = 0;
                _aaid = 0;
            }
            else {
                aaid = 0;
                _aaid = 0;
                baseid = 0;
            }
        }
    },

    /**
    * Set internal flags for all nodes and bonds (internal use)
    * @function setFlag
    */
    setFlag: function (f) {
        for (var i = 0; i < this.atoms.length; ++i)
            this.atoms[i].f = f;
        for (var i = 0; i < this.bonds.length; ++i)
            this.bonds[i].f = f;
    },

    /**
    * Test if the chain contains a node (internal use)
    * @function containsAtom
    */
    containsAtom: function (a) {
        return scil.Utils.indexOf(this.atoms, a) != -1;
    },

    hasBase: function () {
        for (var i = 0; i < this.bases.length; ++i) {
            if (this.bases[i] != null)
                return true;
        }
        return false;
    },

    /**
    * Layout the chain as a line (internal use)
    * @function layoutLine
    */
    layoutLine: function (bondlength) {
        var rect = this.getRect();

        var delta = org.helm.webeditor.bondscale * bondlength;
        var a = this.atoms[0];
        a.p = org.helm.webeditor.Interface.createPoint(rect.left, rect.top);
        for (var i = 1; i < this.atoms.length; ++i) {
            var p = a.p;
            a = this.atoms[i];
            a.p = org.helm.webeditor.Interface.createPoint(p.x + delta, p.y);
        }
    },

    layoutRows: function (bondlength, countperrow) {
        var n = 0;
        var delta = org.helm.webeditor.bondscale * bondlength;
        var x0 = this.atoms[0].p.x;
        for (var i = 1; i < this.atoms.length; ++i) {
            ++n;
            var p = this.atoms[i - 1].p;
            if (n >= countperrow) {
                n = 0;
                //delta = -delta;
                this.bonds[i - 1].z = true;
                this.atoms[i].p = org.helm.webeditor.Interface.createPoint(x0, p.y + Math.abs(delta));
            }
            else {
                this.atoms[i].p = org.helm.webeditor.Interface.createPoint(p.x + delta, p.y);
            }
        }
    },

    /**
    * Layout the chain as a circle (internal use)
    * @function layoutCircle
    */
    layoutCircle: function (bondlength) {
        org.helm.webeditor.Layout.layoutCircle(this.atoms, bondlength, 0);
        //var delta = org.helm.webeditor.bondscale * bondlength;
        //var deg = 360 / (this.atoms.length - 1);
        //var radius = (delta / 2) / Math.sin((deg / 2) * Math.PI / 180);

        //var a = this.atoms[0];
        //a.p = org.helm.webeditor.Interface.createPoint(origin.x + radius, origin.y);
        //for (var i = 1; i < this.atoms.length - 1; ++i)
        //    this.atoms[i].p = this.atoms[i - 1].p.clone().rotateAround(origin, -deg);
    },

    /**
    * Rotate chain (internal use)
    * @function rotate
    */
    rotate: function (deg, origin) {
        if (deg == 0)
            return;

        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            if (origin != null)
                this.atoms[i].p.rotateAround(origin, deg);
            else
                this.atoms[i].p.rotate(deg);

            var a = this.bases[i];
            if (a != null) {
                if (origin != null)
                    a.p.rotateAround(origin, deg);
                else
                    a.p.rotate(deg);
            }
        }
    },

    /**
    * Move chain's coordinates (internal use)
    * @function move
    */
    move: function (delta) {
        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            this.atoms[i].p.offset2(delta);
            var a = this.bases[i];
            if (a != null)
                a.p.offset2(delta);
        }
    },

    /**
    * Layout bases in a RNA chain (internal use)
    * @function layoutBases
    */
    layoutBases: function () {
        var circle = this.isCircle();
        var n = circle ? this.atoms.length - 1 : this.atoms.length;
        for (var i = 0; i < n; ++i) {
            var a = this.bases[i];
            if (a == null)
                continue;

            var center = this.atoms[i];
            var b1 = null;
            var b2 = null;
            if (i == 0) {
                if (circle)
                    b1 = this.bonds[this.bonds.length - 1];
                b2 = this.bonds[i];
            }
            else {
                b1 = this.bonds[i - 1];
                b2 = this.bonds[i];
            }

            if (b1 != null && b2 != null) {
                var a1 = b1.a1 == center ? b1.a2 : b1.a1;
                var a2 = b2.a1 == center ? b2.a2 : b2.a1;

                var ang = center.p.angleAsOrigin(a1.p, a2.p);
                if (Math.abs(ang - 180) > 10)
                    a.p = a1.p.clone().rotateAround(center.p, 180 + ang / 2);
                else
                    a.p = a1.p.clone().rotateAround(center.p, -90);
            }
            else if (b1 != null) {
                var a1 = b1.a1 == center ? b1.a2 : b1.a1;
                a.p = a1.p.clone().rotateAround(center.p, -90);
            }
            else if (b2 != null) {
                var a2 = b2.a1 == center ? b2.a2 : b2.a1;
                a.p = a2.p.clone().rotateAround(center.p, 90);
            }
        }
    },

    /**
    * Get the rectangle of the chain (internal use)
    * @function getRect
    */
    getRect: function () {
        return org.helm.webeditor.Layout.getRect(this.atoms);
    },

    /**
    * Get the sequence of natural analog of the chain
    * @function getSequence
    */
    getSequence: function (highlightselection) {
        var s = "";
        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        var lastbt = null;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            var bt = a.biotype();
            if (bt == org.helm.webeditor.HELM.AA) {
                var mon = org.helm.webeditor.Monomers.getMonomer(a);
                if (mon != null) {
                    if (lastbt != org.helm.webeditor.HELM.AA) {
                        if (s.length > 0 && s.substr(s.length - 1) != "-")
                            s += "-";
                    }
                    var c = scil.Utils.isNullOrEmpty(mon.na) ? "?" : mon.na;
                    if (highlightselection && a.selected)
                        c = "<span style='background:#bbf;'>" + c + "</span>";
                    s += c;
                }
            }
            else if (bt == org.helm.webeditor.HELM.SUGAR) {
                var b = this.bases[i];
                var mon = org.helm.webeditor.Monomers.getMonomer(b);
                if (mon != null) {
                    if (lastbt != org.helm.webeditor.HELM.SUGAR && lastbt != org.helm.webeditor.HELM.LINKER) {
                        if (s.length > 0 && s.substr(s.length - 1) != "-")
                            s += "-";
                    }
                    var c = scil.Utils.isNullOrEmpty(mon.na) ? "?" : mon.na;
                    if (highlightselection && b.selected)
                        c = "<span style='background:#bbf;'>" + c + "</span>";
                    s += c;
                }
            }
            lastbt = bt;
        }

        return s;
    },

    /**
    * Get HELM notation (internal use)
    * @function getHelm
    */
    getHelm: function (ret, highlightselection, m, groupatom) {
        var sequence = [];
        var aaid = 0;
        var firstseqid = null;
        var lastseqid = null;
        var n = this.isCircle() ? this.atoms.length - 1 : this.atoms.length;
        var chn = [];
        var lastbt = null;

        var count = 0;
        for (var i = 0; i < n; ++i) {
            var a = this.atoms[i];
            //if (a.hidden)
            //    continue;

            var bt = a.biotype();
            if (bt == org.helm.webeditor.HELM.LINKER || bt == org.helm.webeditor.HELM.SUGAR || bt == org.helm.webeditor.HELM.NUCLEOTIDE)
                bt = org.helm.webeditor.HELM.BASE;
            if (bt != lastbt || bt == org.helm.webeditor.HELM.CHEM || bt == org.helm.webeditor.HELM.BLOB) {
                var prefix = null;
                if (bt == org.helm.webeditor.HELM.BASE)
                    prefix = "RNA";
                else if (bt == org.helm.webeditor.HELM.AA)
                    prefix = "PEPTIDE";
                else if (bt == org.helm.webeditor.HELM.CHEM)
                    prefix = "CHEM";
                else if (bt == org.helm.webeditor.HELM.BLOB)
                    prefix = a.elem == "Group" ? "G" : "BLOB";

                var seqid = prefix + (++ret.chainid[prefix]);
                if (prefix == "G") {
                    org.helm.webeditor.IO.getGroupHelm(ret, seqid, a, highlightselection);
                    if (this.atoms.length == 1)
                        ret.groupatoms[seqid] = a;
                }

                if (i == 0) {
                    if (a.biotype() == org.helm.webeditor.HELM.SUGAR && a.bio.annotation == "5'ss" || a.bio.annotation == "5'as")
                        ret.annotations[seqid] = { strandtype: a.bio.annotation.substr(2) };
                }

                if (i > 0 && lastseqid != null) {
                    var lasta = this.atoms[i - 1];
                    var b = this.bonds[i - 1];
                    var r1 = b.a1 == lasta ? b.r1 : b.r2;
                    var r2 = b.a1 == lasta ? b.r2 : b.r1;
                    var ratio1 = b.a1 == lasta ? b.ratio1 : b.ratio2;
                    var ratio2 = b.a1 == lasta ? b.ratio2 : b.ratio1;

                    a._aaid = 1;
                    org.helm.webeditor.IO.addConnection(ret, lastseqid, seqid, lasta, a, r1, r2, ratio1, ratio2, null, highlightselection && b.selected);

                    if (!scil.Utils.startswith(lastseqid, "G")) {
                        ++count;
                        ret.sequences[lastseqid] = this._getSequence(sequence, m);
                    }
                    ret.chains[lastseqid] = chn;
                }

                if (firstseqid == null)
                    firstseqid = seqid;

                chn = [];
                aaid = 0;
                sequence = [];
                lastseqid = seqid;
                lastbt = bt;
            }

            a._aaid = ++aaid;
            chn.push(a);

            //if (sequence.length == 0) || aaid > 1 && !(i > 0 && a.biotype() == org.helm.webeditor.HELM.LINKER && this.atoms[i - 1].biotype() == org.helm.webeditor.HELM.SUGAR)) {
            sequence.push({ atoms: [], code: "", type: a.biotype() });
            //}

            var code = org.helm.webeditor.IO.getCode(a, highlightselection);
            sequence[sequence.length - 1].atoms.push(a);

            if (this.bases[i] != null) {
                var b = this.bases[i];
                code += org.helm.webeditor.IO.getCode(b, highlightselection, true);
                sequence[sequence.length - 1].atoms.push(b);
                b._aaid = ++aaid;
                chn.push(b);
            }

            sequence[sequence.length - 1].code += code;
        }

        if (sequence.length > 0) {
            if (!scil.Utils.startswith(lastseqid, "G")) {
                ret.sequences[lastseqid] = this._getSequence(sequence, m);
                if (count == 0 && groupatom != null && !scil.Utils.isNullOrEmpty(groupatom.ratio))
                    ret.ratios[lastseqid] = groupatom.ratio;
            }
            ret.chains[lastseqid] = chn;

            if (groupatom != null && firstseqid == lastseqid && !scil.Utils.isNullOrEmpty(groupatom.tag))
                chn.annotation = groupatom.tag;
        }

        if (this.isCircle()) {
            var b = this.bonds[this.bonds.length - 1];
            //// RNA1,RNA1,1:R1-21:R2
            //var conn;
            //if (this.atoms[0] == b.a1)
            //    conn = b.a1._aaid + ":R" + b.r1 + "-" + b.a2._aaid + ":R" + b.r2;
            //else
            //    conn = b.a2._aaid + ":R" + b.r2 + "-" + b.a1._aaid + ":R" + b.r1;

            //var tag = "";
            //if (!scil.Utils.isNullOrEmpty(b.tag))
            //    tag = '\"' + b.tag.replace(/"/g, "\\\"") + '\"';

            //ret.connections.push(firstseqid + "," + lastseqid + "," + conn);

            org.helm.webeditor.IO.addConnection(ret, firstseqid, lastseqid, b.a1, b.a2, b.r1, b.r2, null, null, b.tag, highlightselection && b.selected);
        }
    },

    _getSequence: function (sequence, m) {
        var ret = "";
        var needseparator = false;
        for (var i = 0; i < sequence.length; ++i) {
            if (i > 0) {
                if (needseparator || !(sequence[i - 1].type == org.helm.webeditor.HELM.SUGAR && sequence[i].type == org.helm.webeditor.HELM.LINKER)) {
                    ret += ".";
                    needseparator = false;
                }
            }

            var br = this._getBracket(sequence[i].atoms[0], m);
            var repeat = br == null ? null : br.getSubscript(m);

            if (!scil.Utils.isNullOrEmpty(repeat)) {
                var end = -1;
                var atoms = scil.clone(br.atoms);
                for (var k = i; k < sequence.length; ++k) {
                    var list = sequence[k].atoms;
                    for (var j = 0; j < list.length; ++j) {
                        var p = scil.Utils.indexOf(atoms, list[j]);
                        if (p >= 0) {
                            atoms.splice(p, 1);
                        }
                        else {
                            // ERROR: invalid bracket
                            end = -2;
                            break;
                        }
                    }

                    if (end == -2)
                        break;
                    if (atoms.length == 0) {
                        end = k;
                        break;
                    }
                }

                if (end >= 0) {
                    if (ret.length > 0 && !scil.Utils.endswith(ret, "."))
                        ret += ".";

                    if (i == end) {
                        ret += sequence[i].code;
                    }
                    else {
                        ret += "(";
                        var st = i;
                        for (var k = i; k <= end; ++k) {
                            if (k > st) {
                                if (!(sequence[k - 1].type == org.helm.webeditor.HELM.SUGAR && sequence[k].type == org.helm.webeditor.HELM.LINKER))
                                    ret += ".";
                            }
                            ret += sequence[k].code;
                        }
                        ret += ")";
                    }
                    ret += "'" + repeat + "'";

                    i = end;
                    needseparator = true;
                    continue;
                }
            }

            ret += sequence[i].code;
        }
        return ret;
    },

    _getBracket: function (a, m) {
        if (m == null)
            return null;

        for (var k = 0; k < m.graphics.length; ++k) {
            var br = JSDraw2.Bracket.cast(m.graphics[k]);
            if (br == null)
                continue;

            if (scil.Utils.indexOf(br.atoms, a) >= 0)
                return br;
        }

        return null;
    },

    /**
    * Get a node by its Monomer ID (internal use)
    * @function getAtomByAAID
    */
    getAtomByAAID: function (aaid) {
        if (!(aaid > 0))
            return null;

        for (var i = 0; i < this.atoms.length; ++i) {
            if (this.atoms[i]._aaid == aaid)
                return this.atoms[i];
        }
        for (var i = 0; i < this.bases.length; ++i) {
            if (this.bases[i]._aaid == aaid)
                return this.bases[i];
        }

        return null;
    }
});

scil.apply(org.helm.webeditor.Chain, {
    /**
    * Get the chain of a given node (internal use)
    * @function getChain
    */
    getChain: function (m, startatom) {
        if (startatom == null)
            return null;
        var chains = this._getChains(m, startatom);
        return chains == null ? null : chains[0];
    },

    /**
    * Get all chains (internal use)
    * @function getChains
    */
    getChains: function (m, branchcollection) {
        return this._getChains(m, null, branchcollection);
    },

    /**
    * Set Chain ID (internal use)
    * @function _setChainID
    */
    _setChainID: function (chain, chainid) {
        for (var k = 0; k < chain.atoms.length; ++k) {
            chain.atoms[k]._chainid = chainid;
            if (chain.bases[k] != null)
                chain.bases[k]._chainid = chainid;
        }
    },

    /**
    * Remove Chain ID (internal use)
    * @function _removeChainID
    */
    _removeChainID: function (atoms) {
        for (var i = 0; i < atoms.length; ++i)
            delete atoms[i]._chainid;
    },

    /**
    * Inner loop getting chains (internal use)
    * @function _getChains
    */
    _getChains: function (m, startatom, branchcollection) {
        var b0 = null;
        var bonds = [];
        var branches = [];
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.r1 == 1 && b.r2 == 2 || b.r1 == 2 && b.r2 == 1) {
                bonds.push(b);
                if (b0 == null && startatom != null && b.a1 == startatom || b.a2 == startatom)
                    b0 = bonds.length - 1;
            }
            else {
                branches.push(b);
            }
        }

        if (startatom != null && b0 == null)
            return null;

        var chains = [];
        while (bonds.length > 0) {
            var chain = new org.helm.webeditor.Chain();
            chains.splice(0, 0, chain);

            var b = null;
            if (b0 == null) {
                b = bonds[bonds.length - 1];
                bonds.splice(bonds.length - 1, 1);
            }
            else {
                b = bonds[b0];
                bonds.splice(b0, 1);
                b0 = null;
            }
            
            var head = b.r1 == 2 ? b.a1 : b.a2;
            var tail = b.r2 == 1 ? b.a2 : b.a1;

            chain.bonds.push(b);
            chain.atoms.push(head);
            chain.atoms.push(tail);

            while (bonds.length > 0) {
                var found = 0;
                for (var i = bonds.length - 1; i >= 0; --i) {
                    b = bonds[i];
                    if (b.a1 == head || b.a2 == head) {
                        bonds.splice(i, 1);
                        head = b.a1 == head ? b.a2 : b.a1;
                        chain.bonds.splice(0, 0, b);
                        chain.atoms.splice(0, 0, head);

                        ++found;
                    }
                    else if (b.a1 == tail || b.a2 == tail) {
                        bonds.splice(i, 1);
                        tail = b.a1 == tail ? b.a2 : b.a1;
                        chain.bonds.push(b);
                        chain.atoms.push(tail);

                        ++found;
                    }
                }

                if (found == 0)
                    break;
            }

            if (startatom != null)
                break;
        }

        m.clearFlag();
        for (var i = 0; i < chains.length; ++i) {
            var atoms = chains[i].atoms;
            for (var k = 0; k < atoms.length; ++k)
                atoms[k].f = true;
        }

        if (startatom == null) {
            var atoms = m.atoms;
            for (var k = 0; k < atoms.length; ++k) {
                var a = atoms[k];
                if (a.f)
                    continue;

                //if (a.biotype() == org.helm.webeditor.HELM.AA || a.biotype() == org.helm.webeditor.HELM.SUGAR || a.biotype() == org.helm.webeditor.HELM.LINKER) {
                if (a.biotype() != org.helm.webeditor.HELM.BASE) {
                    a.f = true;
                    var chain = new org.helm.webeditor.Chain();
                    chains.splice(0, 0, chain);
                    chain.atoms.push(a);
                }
            }
        }

        for (var i = 0; i < chains.length; ++i) {
            var atoms = chains[i].atoms;
            var bonds = chains[i].bonds;

            if (chains[i].isCircle()) {
                // rotate circle if the first atom is a linker (P)
                if (atoms[0].biotype() == org.helm.webeditor.HELM.LINKER && atoms[1].biotype() == org.helm.webeditor.HELM.SUGAR) {
                    atoms.splice(0, 1);
                    atoms.push(atoms[0]);

                    bonds.push(bonds[0]);
                    bonds.splice(0, 1);
                }

                // rotate if RNA/PEPTIDE/CHEM circle
                for (var j = 0; j < atoms.length - 1; ++j) {
                    var bt1 = atoms[j].biotype();
                    if (bt1 == org.helm.webeditor.HELM.LINKER)
                        bt1 = org.helm.webeditor.HELM.SUGAR;

                    var bt2 = atoms[j + 1].biotype();
                    if (bt2 == org.helm.webeditor.HELM.LINKER)
                        bt2 = org.helm.webeditor.HELM.SUGAR;

                    if (bt1 != bt2) {
                        for (var k = 0; k <= j; ++k) {
                            atoms.splice(0, 1);
                            atoms.push(atoms[0]);

                            bonds.push(bonds[0]);
                            bonds.splice(0, 1);
                        }
                        break;
                    }
                }
            }

            // detect bases
            for (var k = 0; k < atoms.length; ++k) {
                var a = atoms[k];
                if (a.biotype() == org.helm.webeditor.HELM.SUGAR) {
                    for (var j = branches.length - 1; j >= 0; --j) {
                        var at = null;
                        var b = branches[j];
                        if (b.a1 == a && b.r1 == 3 && b.r2 == 1 && b.a2.biotype() == org.helm.webeditor.HELM.BASE)
                            at = chains[i].bases[k] = b.a2;
                        else if (b.a2 == a && b.r2 == 3 && b.r1 == 1 && b.a1.biotype() == org.helm.webeditor.HELM.BASE)
                            at = chains[i].bases[k] = b.a1;

                        if (at != null) {
                            chains[i].basebonds.push(b);
                            branches.splice(j, 1);
                            at.f = true;
                        }
                    }
                }
            }
        }

        if (branchcollection != null) {
            var list = [];
            for (var i = 0; i < branches.length; ++i) {
                var b = branches[i];
                if (!b.a1.f) {
                    b.a1.f = true;
                    list.push(b.a1);
                }
                if (!b.a2.f) {
                    b.a2.f = true;
                    list.push(b.a2);
                }
            }

            branchcollection.bonds = branches;
            branchcollection.atoms = list;
        }

        return chains;
    }
});﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* Layout class
* @class org.helm.webeditor.Layout
*/
org.helm.webeditor.Layout = {
    /**
    * Clean/lay out a molecule
    * @function clean
    */
    clean: function (jsd, a) {
        var m = jsd.m;
        var bondlength = jsd.bondlength;

        m.clearFlag();
        var chains = org.helm.webeditor.Chain._getChains(m, a);

        org.helm.webeditor.Chain._removeChainID(m.atoms);
        for (var i = 0; i < chains.length; ++i) {
            var chain = chains[i];

            // set chain id
            org.helm.webeditor.Chain._setChainID(chain, i);

            if (!chain.hasBase() && chain.atoms.length > 50) {
                chain.layoutRows(bondlength, 10);
            }
            else {
                if (chain.isCircle()) {
                    chain.layoutCircle(bondlength);
                }
                else {
                    if (!this.layoutInnerCircle(chain, bondlength, m, i))
                        chain.layoutLine(bondlength);
                }
                chain.layoutBases();
            }

            chain.setFlag(true);
            chain.resetIDs();
        }

        jsd.updateGroupRect();

        this.layoutCrossChainBonds(m, chains, bondlength);
        this.layoutBranches(m);

        this.layoutFragments(m, bondlength);

        // clear chain id
        org.helm.webeditor.Chain._removeChainID(m.atoms);
    },

    layoutFragments: function (m, bondlength) {
        var frags = m.splitFragments(true);
        if (frags.length < 2)
            return;

        var r0 = frags[0].rect();
        var x = r0.center().x;
        var y = r0.bottom() + 3 * bondlength;
        for (var i = 1; i < frags.length; ++i) {
            var frag = frags[i];
            var r = frag.rect();
            frag.offset(x - r.center().x, y - r.top);
            y += r.height + 3 * bondlength;
        }
    },

    /**
    * Reset Monomer IDs (internal use)
    * @function resetIDs
    */
    resetIDs: function (m, aaid) {
        var chains = org.helm.webeditor.Chain._getChains(m);
        for (var i = 0; i < chains.length; ++i)
            chains[i].resetIDs(aaid);
    },

    /**
    * Lay out inner circle (internal use)
    * @function layoutInnerCircle
    */
    layoutInnerCircle: function (chain, bondlength, m, chainid) {
        var pairs = [];
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.a1._chainid != null && b.a2._chainid != null && b.a1._chainid == chainid && b.a2._chainid == chainid && scil.Utils.indexOf(chain.bonds, b) < 0 && scil.Utils.indexOf(chain.basebonds, b) < 0) {
                var ai1 = scil.Utils.indexOf(chain.atoms, b.a1);
                var ai2 = scil.Utils.indexOf(chain.atoms, b.a2);
                var p1 = { a1: ai1 < ai2 ? ai1 : ai2, a2: ai1 < ai2 ? ai2 : ai1 };
                pairs.push(p1);
            }
        }

        if (pairs.length == 0)
            return false;

        // find the biggest circle
        var pair = pairs[0];
        for (var i = 1; i < pairs.length; ++i) {
            var r = pairs[i];
            if (r.a1 >= pair.a1 && r.a1 <= pair.a2 && r.a1 >= pair.a1 && r.a1 <= pair.a2) {
                pair = r;
            }
            else if (pair.a1 < r.a1 || pair.a1 > r.a2 || pair.a2 < r.a1 || pair.a2 > r.a2) {
                if (r.a2 - r.a1 > pair.a2 - pair.a1)
                    pair = r;
            }
        }

        var atoms = [];
        for (var i = pair.a1; i <= pair.a2; ++i)
            atoms.push(chain.atoms[i]);
        atoms.push(atoms[0]);
        this.layoutCircle(atoms, bondlength, -360 / (atoms.length - 1) / 2);

        var delta = org.helm.webeditor.bondscale * bondlength;
        var p = chain.atoms[pair.a1].p.clone();
        for (var i = pair.a1 - 1; i >= 0; --i) {
            p.x += delta;
            chain.atoms[i].p = p.clone();
        }

        p = chain.atoms[pair.a2].p.clone();
        for (var i = pair.a2 + 1; i < chain.atoms.length; ++i) {
            p.x += delta;
            chain.atoms[i].p = p.clone();
        }

        return true;
    },

    /**
    * Lay out circle (internal use)
    * @function layoutCircle
    */
    layoutCircle: function (atoms, bondlength, startdeg) {
        var rect = this.getRect(atoms);
        var origin = rect.center();

        var delta = org.helm.webeditor.bondscale * bondlength;
        var deg = 360 / (atoms.length - 1);
        var radius = (delta / 2) / Math.sin((deg / 2) * Math.PI / 180);

        var a = atoms[0];
        a.p = org.helm.webeditor.Interface.createPoint(origin.x + radius, origin.y);
        if (startdeg != null && startdeg != 0)
            a.p.rotateAround(origin, startdeg);

        for (var i = 1; i < atoms.length - 1; ++i)
            atoms[i].p = atoms[i - 1].p.clone().rotateAround(origin, -deg);
    },

    /**
    * Lay out cross-chain bonds (internal use)
    * @function layoutCrossChainBonds
    */
    layoutCrossChainBonds: function (m, chains, bondlength) {
        var fixed = {};
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.a1._chainid != null && b.a2._chainid != null && b.a1._chainid != b.a2._chainid) {
                var a1, a2;
                if (fixed[b.a1._chainid] && fixed[b.a2._chainid]) {
                    continue;
                }
                else if (fixed[b.a1._chainid]) {
                    a1 = b.a1;
                    a2 = b.a2;
                }
                else if (fixed[b.a2._chainid]) {
                    a1 = b.a2;
                    a2 = b.a1;
                }
                else {
                    var chain1 = chains[b.a1._chainid];
                    var chain2 = chains[b.a2._chainid];
                    if (chain1.atoms.length < chain2.atoms.length) {
                        a1 = b.a2;
                        a2 = b.a1;
                    }
                    else {
                        a1 = b.a1;
                        a2 = b.a2;
                    }
                }

                if (b.type == JSDraw2.BONDTYPES.UNKNOWN) {
                    // hydrogen bond
                    if (b.a1.p.y > b.a2.p.y) {
                        a2 = b.a1;
                        a1 = b.a2;
                    }
                    else {
                        a2 = b.a2;
                        a1 = b.a1;
                    }
                    var chain = chains[a2._chainid];
                    chain.rotate(180);

                    var delta = a1.p.clone().offset(0, bondlength * org.helm.webeditor.bondscale).offset(-a2.p.x, -a2.p.y);
                    chain.move(delta);
                }
                else {
                    // cross-chain connection
                    var bonds = m.getNeighborBonds(a1);
                    if (bonds.length == 3) {
                        scil.Utils.delFromArray(bonds, b);

                        var p1 = bonds[0].otherAtom(a1).p;
                        var p2 = bonds[1].otherAtom(a1).p;
                        var p = new JSDraw2.Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
                        if (p.distTo(a1.p) < bondlength / 30) {
                            // p1, a1.p and p2 in a line
                            p = p1.clone();
                            p.rotateAround(a1.p, -90, bondlength * 3);
                        }
                        else {
                            p.rotateAround(a1.p, 180, bondlength * 3);
                        }

                        var chain = chains[a2._chainid];
                        p.offset(-a2.p.x, -a2.p.y);
                        chain.move(p);

                        bonds = m.getNeighborBonds(a2);
                        if (bonds.length == 3) {
                            scil.Utils.delFromArray(bonds, b);

                            var deg;
                            var c = a2.p.clone();

                            var ang1 = a1.p.angleTo(c);

                            var p1 = bonds[0].otherAtom(a2).p;
                            var p2 = bonds[1].otherAtom(a2).p;
                            var p = new JSDraw2.Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
                            if (p.distTo(c) < bondlength / 30) {
                                // p1, a2.p and p2 in a line
                                var ang2 = p2.angleTo(c);
                                deg = (ang1 - ang2) - 90;
                            }
                            else {
                                var ang2 = p.angleTo(c);
                                deg = (ang1 + 180) - ang2;
                            }

                            chain.rotate(deg, c);
                        }
                    }
                    else {
                        var p = a1.p.clone().offset(0, bondlength * 3);
                        if (this._hasOverlap(chains, fixed, p, bondlength / 30))
                            p = a1.p.clone().offset(0, -bondlength * 3);
                        var delta = p.offset(-a2.p.x, -a2.p.y);
                        chains[a2._chainid].move(delta);
                    }
                }

                fixed[a1._chainid] = true;
                fixed[a2._chainid] = true;
            }
        }
    },

    _hasOverlap: function (chains, fixed, p, tor) {
        for (var i = 0; i < chains.length; ++i) {
            if (!fixed[i])
                continue;

            var atoms = chains[i].atoms;
            for (var k = 0; k < atoms.length; ++k) {
                if (atoms[k].p.distTo(p) < tor)
                    return true;
            }
        }

        return false;
    },

    /**
    * Layout branches (internal use)
    * @function layoutBranches
    */
    layoutBranches: function (m) {
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (!b.f && b.a1.f != b.a2.f) {
                var center = b.a1.f ? b.a1 : b.a2;
                var a = b.a1.f ? b.a2 : b.a1;

                var b1 = null;
                var b2 = null;
                var bonds = m.getNeighborBonds(center);
                for (var k = bonds.length - 1; k >= 0; --k) {
                    var n = bonds[k];
                    if (n.f) {
                        if (b1 == null && n.a1 == center && n.r1 == 2 || n.a2 == center && n.r2 == 2) {
                            b1 = n;
                            bonds.splice(i, 0);
                        }
                        else if (b2 == null && n.a1 == center && n.r1 == 1 || n.a2 == center && n.r2 == 1) {
                            b2 = n;
                            bonds.splice(i, 0);
                        }
                    }
                }

                if (b1 != null || b2 != null) {
                    if (b1 != null && b2 != null) {
                        var a1 = b1.a1 == center ? b1.a2 : b1.a1;
                        var a2 = b2.a1 == center ? b2.a2 : b2.a1;

                        var ang = center.p.angleAsOrigin(a1.p, a2.p);
                        if (Math.abs(ang - 180) > 10)
                            a.p = a1.p.clone().rotateAround(center.p, ang / 2);
                        else
                            a.p = a1.p.clone().rotateAround(center.p, 90);
                    }
                    else {
                        if (b1 != null) {
                            var a1 = b1.a1 == center ? b1.a2 : b1.a1;
                            a.p = a1.p.clone().rotateAround(center.p, 180);
                        }
                        else if (b2 != null) {
                            var a2 = b2.a1 == center ? b2.a2 : b2.a1;
                            a.p = a2.p.clone().rotateAround(center.p, 180);
                        }
                    }

                    b.f = b.a1.f = b.a2.f = true;
                }
            }
        }
    },

    /**
    * Get rectangle (internal use)
    * @function getRect
    */
    getRect: function (atoms) {
        var a = atoms[0];
        var x1 = a.p.x;
        var y1 = a.p.y;
        var x2 = x1;
        var y2 = y1;

        for (var i = 1; i < atoms.length; ++i) {
            var p = atoms[i].p;
            if (p.x < x1)
                x1 = p.x;
            else if (p.x > x2)
                x2 = p.x;
            if (p.y < y1)
                y1 = p.y;
            else if (p.y > y2)
                y2 = p.y;
        }

        return org.helm.webeditor.Interface.createRect(x1, y1, x2 - x1, y2 - y1);
    }
};
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* IO class
* @class org.helm.webeditor.IO
*/
org.helm.webeditor.IO = {
    kVersion: "V2.0",

    /**
    * Get HELM Notation
    * @function getHelm
    */
    getHelm: function (m, highlightselection) {
        // I#12164
        for (var i = 0; i < m.atoms.length; ++i) {
            var a = m.atoms[i];
            if (!org.helm.webeditor.isHelmNode(a) && !a.hidden)
                return null;
        }

        var ret = { chainid: { RNA: 0, PEPTIDE: 0, CHEM: 0, BLOB: 0, G: 0 }, sequences: {}, connections: [], chains: {}, pairs: [], groupatoms: [], groups: {}, annotations: {}, singletons: {}, ratios: {} };
        this.getHelm2(m, highlightselection, ret);

        for (var k in ret.chains) {
            var chain = ret.chains[k];
            var a = chain[0];
            if (chain.length == 1 && a.biotype() == org.helm.webeditor.HELM.BLOB && !scil.Utils.isNullOrEmpty(a.tag))
                chain.annotation = a.tag;
        }

        for (var k in ret.groups) {
            if (ret.groups[k].length == 1)
                ret.singletons[k] = ret.groups[k][0];
        }

        for (var k in ret.singletons)
            delete ret.groups[k];

        return this.getHelmString(ret, highlightselection);
    },

    getHelm2: function (m, highlightselection, ret, groupatom) {
        var branches = {};
        var chains = org.helm.webeditor.Chain.getChains(m, branches);

        for (var i = 0; i < m.atoms.length; ++i)
            m.atoms[i]._aaid = null;

        for (var i = 0; i < chains.length; ++i) {
            var chain = chains[i];
            chain.getHelm(ret, highlightselection, m, groupatom);
        }

        for (var i = 0; i < branches.atoms.length; ++i) {
            var a = branches.atoms[i];
            if (a.biotype() == org.helm.webeditor.HELM.CHEM) {
                var id = "CHEM" + (++ret.chainid.CHEM);
                ret.sequences[id] = this.getCode(a, highlightselection);
                ret.chains[id] = [a];
                a._aaid = 1;
            }
            else if (a.biotype() == org.helm.webeditor.HELM.BLOB) {
                var id;
                if (a.elem == "Group") {
                    id = "G" + (++ret.chainid.G);
                    this.getGroupHelm(ret, id, a, highlightselection);
                }
                else {
                    id = "BLOB" + (++ret.chainid.BLOB);
                    ret.sequences[id] = this.getCode(a, highlightselection);
                }
                ret.chains[id] = [a];
                a._aaid = 1;
            }
            else {
                // error
                return null;
            }
        }

        var groups = [];
        for (var i = 0; i < m.graphics.length; ++i) {
            var g = JSDraw2.Group.cast(m.graphics[i]);
            if (g != null)
                groups.push(g);
        }

        for (var i = 0; i < groups.length; ++i) {
            var g = groups[i];
            if (scil.Utils.isNullOrEmpty(g.tag))
                continue;

            for (var c in ret.chains) {
                if (this.allBelongToGroup(ret.chains[c], g)) {
                    ret.chains[c].annotation = g.tag;
                    break;
                }
            }
        }

        for (var id in ret.groupatoms) {
            var a = ret.groupatoms[id];
            for (var i = 0; i < groups.length; ++i) {
                var g = groups[i];
                if (g.a == a) {
                    groups.splice(i, 1);
                    this._scanGroup(ret, g, id);
                    break;
                }
            }
        }

        var groupids = [];
        for (var i = 0; i < groups.length; ++i) {
            var g = groups[i];

            var prefix = "G";
            var id = prefix + (++ret.chainid[prefix]);
            groupids[i] = id;

            this._scanGroup(ret, g, id);
        }

        for (var i = 0; i < groups.length; ++i) {
            var g = groups[i];
            var combo = ret.groups[groupids[i]];
            if (combo == null)
                continue;

            for (var j = 0; j < groups.length; ++j) {
                var g2 = groups[j];
                if (g2.group == g) {
                    var id = groupids[j];
                    if (scil.Utils.indexOf(combo, id) < 0)
                        combo.push(id);
                }
            }

            for (var id in ret.groupatoms) {
                var a = ret.groupatoms[id];
                if (a.group == g) {
                    if (scil.Utils.indexOf(combo, id) < 0)
                        combo.push(id);
                    break;
                }
            }
        }

        // RNA1,RNA2,5:pair-11:pair
        for (var i = 0; i < branches.bonds.length; ++i) {
            var b = branches.bonds[i];

            var tag = "";
            if (!scil.Utils.isNullOrEmpty(b.tag))
                tag = '\"' + b.tag.replace(/"/g, "\\\"") + '\"';

            if (b.type == JSDraw2.BONDTYPES.UNKNOWN) {
                var c1 = this.findChainID(ret.chains, b.a1);
                var c2 = this.findChainID(ret.chains, b.a2);
                var s = c1 + "," + c2 + "," + b.a1._aaid + ":pair-" + b.a2._aaid + ":pair";
                ret.pairs.push(s + tag);
            }
            else {
                var c1 = this.findChainID(ret.chains, b.a1);
                var c2 = this.findChainID(ret.chains, b.a2);
                this.addConnection(ret, c1, c2, b.a1, b.a2, b.r1, b.r2, b.ratio1, b.ratio2, b.tag, highlightselection && b.selected);
            }
        }
    },

    _scanGroup: function (ret, g, id) {
        var combo = [];
        for (var c in ret.chains) {
            if (this.allBelongToGroup(ret.chains[c], g)) {
                combo.push(c);
                if (!scil.Utils.isNullOrEmpty(g.ratio))
                    ret.ratios[c] = g.ratio;
            }
        }

        ret.groups[id] = combo;
        ret.ratios[id] = g.ratio;
    },

    getGroupHelm: function (ret, id, a, highlightselection) {
        var existing = {};
        for (var k in ret.sequences)
            existing[k] = true;

        var m = a.superatom;
        if (a.hidden) {
            ret.groupatoms[id] = a;
        }
        else {
            this.getHelm2(m, highlightselection, ret, a);

            var combo = [];
            for (var k in ret.sequences) {
                if (!existing[k])
                    combo.push(k); //{ chain: k, ratio: a.ratio });
            }

            ret.groups[id] = combo;
        }
    },

    addConnection: function (ret, c1, c2, a1, a2, r1, r2, ratio1, ratio2, tag, h) {
        var ai1 = a1.biotype() == org.helm.webeditor.HELM.BLOB ? "?" : a1._aaid;
        var ai2 = a2.biotype() == org.helm.webeditor.HELM.BLOB ? "?" : a2._aaid;
        ret.connections.push({ c1: c1, c2: c2, ai1: ai1, ai2: ai2, r1: r1, r2: r2, ratio1: ratio1, ratio2: ratio2, tag: tag, h: h });
    },

    renderConnection: function (ret, conn) {
        // if it is G1(PEPTID1), then directly use PEPTIDE1, and not G1
        if (ret.singletons[conn.c1] != null) {
            conn.c1 = ret.singletons[conn.c1];
            if (conn.ai1 > 0)
                conn.ai1 = "?";
        }
        if (ret.singletons[conn.c2] != null) {
            conn.c2 = ret.singletons[conn.c2];
            if (conn.ai2 > 0)
                conn.ai2 = "?";
        }

        var c = conn.c1 + "," + conn.c2;
        c += "," + org.helm.webeditor.IO.connectionStr(conn.ai1, conn.r1, conn.ai2, conn.r2);

        if (!scil.Utils.isNullOrEmpty(conn.tag))
            c += '\"' + conn.tag.replace(/"/g, "\\\"") + '\"';

        if (conn.h)
            c = "<span style='background:#bbf;'>" + c + "</span>";
        return c;
    },

    connectionStr: function (aaid1, r1, aaid2, r2) {
        return this.rStr(aaid1, r1) + "-" + this.rStr(aaid2, r2);
    },

    rStr: function (aaid, r) {
        if (typeof (r) == "string" && r.indexOf(':') > 0)
            return r;

        var s = aaid + ":";
        if (r == "*" || r == "?")
            return s + r;
        return s + "R" + r;
    },

    allBelongToGroup: function (atoms, g) {
        for (var i = 0; i < atoms.length; ++i) {
            if (atoms[i].group != g)
                return false;
        }
        return true;
    },

    getHelmString: function (ret, highlightselection) {
        var s = "";
        var keys = [];
        for (var k in ret.sequences)
            keys.push(k);
        keys.sort();
        for (var i = 0; i < keys.length; ++i) {
            var k = keys[i];
            s += (s == "" ? "" : "|") + k + "{" + ret.sequences[k] + "}";
            var chain = ret.chains[k];
            if (chain != null && !scil.Utils.isNullOrEmpty(chain.annotation))
                s += this.wrapAnnotation(chain.annotation);
        }

        if (s == "")
            return s;

        var count = 0;

        s += "$";
        var groups = [];
        for (var i = 0; i < ret.connections.length; ++i) {
            var c = ret.connections[i];
            s += (++count > 1 ? "|" : "") + this.renderConnection(ret, c);
            if ((c.ratio1 > 0 || c.ratio1 == "?") && (c.ratio2 > 0 || c.ratio1 == "?")) {
                var s2 = c.c1 + ":" + c.ratio1 + "+" + c.c2 + ":" + c.ratio2;
                id = "G" + (++ret.chainid.G);
                groups.push(id + "(" + s2 + ")");
            }
        }
        for (var i = 0; i < ret.pairs.length; ++i)
            s += (++count > 1 ? "|" : "") + ret.pairs[i];

        s += "$";
        var list = [];
        for (var id in ret.groups) {
            var s2 = "";
            var list = ret.groups[id];
            for (var i = 0; i < list.length; ++i) {
                var c = list[i];
                if (ret.singletons[c] != null)
                    c = ret.singletons[c];
                var ratio = ret.ratios[c];
                var separator = ret.ratios[id] == "or" ? "," : "+";
                s2 += (i > 0 ? separator : "") + c + (scil.Utils.isNullOrEmpty(ratio) ? "" : ":" + ratio);
            }
            groups.push(id + "(" + s2 + ")");
        }
        for (var i = 0; i < groups.length; ++i)
            s += (i > 0 ? "|" : "") + groups[i];

        s += "$";

        //RNA1{R(C)P.R(A)P.R(T)}$$$RNA1{ss}$
        var ann = scil.Utils.json2str(ret.annotations, null, true);
        s += ann == "null" || ann == "{}" ? "" : ann;

        s += "$";
        return s + this.kVersion;
    },

    /**
    * Get the natural sequence of the molecule
    * @function getSequence
    */
    getSequence: function (m, highlightselection) {
        var branches = {};
        var chains = org.helm.webeditor.Chain.getChains(m, branches);
        if (chains == null)
            return null;

        var s = "";
        for (var i = 0; i < chains.length; ++i) {
            var s2 = chains[i].getSequence(highlightselection);
            if (scil.Utils.isNullOrEmpty(s2))
                continue;
            if (s != "")
                s += "\r\n";
            s += s2;
        }

        return s;
    },

    /**
    * Get XHELM
    * @function getXHelm
    */
    getXHelm: function (m) {
        var s = this.getHelm(m);
        if (scil.Utils.isNullOrEmpty(s))
            return s;

        var s = "<Xhelm>\n<HelmNotation>" + scil.Utils.escXmlValue(s) + "</HelmNotation>\n";

        var list = this.getMonomers(m);
        if (list != null) {
            s += "<Monomers>\n";
            for (var i in list)
                s += org.helm.webeditor.Monomers.writeOne(list[i]);
            s += "</Monomers>\n";
        }
        s += "</Xhelm>";
        return s;
    },

    /**
    * Get all monomers of a molecule
    * @function getMonomers
    */
    getMonomers: function (m) {
        var ret = {};
        var atoms = m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            var biotype = a.biotype();
            if (!org.helm.webeditor.isHelmNode(a))
                continue;

            var m = org.helm.webeditor.Monomers.getMonomer(a);
            var t = biotype + " " + a.elem;
            if (ret[t] == null) {
                var type = null;
                var mt = null;
                if (biotype == org.helm.webeditor.HELM.CHEM) {
                    type = "CHEM";
                    mt = "Undefined"
                }
                else if (biotype == org.helm.webeditor.HELM.AA) {
                    type = "PEPTIDE";
                    mt = "Undefined"
                }
                else if (biotype == org.helm.webeditor.HELM.SUGAR) {
                    type = "RNA";
                    mt = "Backbone"
                }
                else if (biotype == org.helm.webeditor.HELM.BASE) {
                    type = "RNA";
                    mt = "Branch"
                }
                else if (biotype == org.helm.webeditor.HELM.LINKER) {
                    type = "RNA";
                    mt = "Backbone"
                }
                ret[t] = { id: a.elem, mt: mt, type: type, m: m };
            }
        }

        return ret;
    },

    /**
    * Get HELM Code of a monomer (internal use)
    * @function getCode
    */
    getCode: function (a, highlightselection, bracket) {
        var s;
        var blob = false;
        if (typeof (a) == "object" && a.biotype() == org.helm.webeditor.HELM.BLOB) {
            blob = true;
            if (a.elem == "Group") {
                s = "Group";
            }
            else {
                s = a.bio != null && scil.Utils.isNullOrEmpty(a.bio.blobtype) ? "" : a.bio.blobtype;
            }
        }
        else {
            if (typeof (a) == "string") {
                s = a;
            }
            else {
                var m = org.helm.webeditor.Monomers.getMonomer(a);
                if (m.issmiles)
                    s = m.smiles;
                else
                    s = a.elem;
            }

            if (s == "?" && a.bio != null)
                s = a.bio.ambiguity == null ? "?" : a.bio.ambiguity;
            else if (s.length > 1)
                s = "[" + s + "]";
        }

        if (!blob)
            s += this.wrapAnnotation(a.tag);

        if (bracket)
            s = "(" + s + ")";
        if (highlightselection && a.selected)
            s = "<span style='background:#bbf;'>" + s + "</span>";
        return s;
    },

    wrapAnnotation: function (s) {
        if (!scil.Utils.isNullOrEmpty(s))
            return '\"' + s.replace(/"/g, "\\\"") + '\"';
        return "";
    },

    /**
    * Find the chain ID based on monomer (internal use)
    * @function findChainID
    */
    findChainID: function (chains, a) {
        for (var k in chains) {
            var atoms = chains[k];
            if (scil.Utils.indexOf(atoms, a) >= 0)
                return k;
        }
        return null;
    },

    /**
    * Read a generic string (internal use)
    * @function read
    */
    read: function (plugin, s, format, renamedmonomers, sugar, linker, separator) {
        if (scil.Utils.isNullOrEmpty(s))
            return 0;

        var s2 = s.toUpperCase();
        if (scil.Utils.isNullOrEmpty(format)) {
            if (/^((RNA)|(PEPTIDE)|(CHEM)|(BLOB))[0-9]+/.test(s2))
                format = "HELM";
            else if (/^[A|G|T|C|U]+[>]?$/.test(s2))
                format = "RNA";
            else if (/^[A|C-I|K-N|P-T|V|W|Y|Z]+[>]?$/.test(s2))
                format = "Peptide";
            else
                throw "Cannot detect the format using nature monomer names";
        }

        var origin = org.helm.webeditor.Interface.createPoint(0, 0);
        if (format == "HELM") {
            return this.parseHelm(plugin, s, origin, renamedmonomers);
        }
        else if (format == "Peptide") {
            var chain = new org.helm.webeditor.Chain();
            var circle = scil.Utils.endswith(s, ">");
            if (circle)
                s = s.substr(0, s.length - 1);
            var ss = this.splitChars(s, separator);
            if (circle)
                ss.push(">");
            return this.addAAs(plugin, ss, chain, origin);
        }
        else if (format == "RNA") {
            var chain = new org.helm.webeditor.Chain();
            var ss = this.splitChars(s, separator);
            return this.addRNAs(plugin, ss, chain, origin, sugar, linker);
        }

        return 0;
    },

    /**
    * Parse a HELM string (internal use)
    * @function parseHelm
    */
    parseHelm: function (plugin, s, origin, renamedmonomers) {
        var n = 0;
        var sections = this.split(s, '$');
        var chains = {};

        var gi = 100;
        var groups = {};
        var groupannotations = {};

        // sequence
        s = sections[0];
        if (!scil.Utils.isNullOrEmpty(s)) {
            var seqs = this.split(s, '|');
            for (var i = 0; i < seqs.length; ++i) {
                var e = this.detachAnnotation(seqs[i]);
                s = e.str;

                p = s.indexOf("{");
                var sid = s.substr(0, p);
                var type = sid.replace(/[0-9]+$/, "").toUpperCase();
                var id = parseInt(sid.substr(type.length));

                var chain = new org.helm.webeditor.Chain(sid);
                chains[sid] = chain;
                chain.type = type;

                s = s.substr(p + 1);
                p = s.indexOf('}');
                s = s.substr(0, p);

                var n2 = 0;
                var ss = this.split(s, '.');
                if (type == "PEPTIDE")
                    n2 = this.addAAs(plugin, ss, chain, origin, renamedmonomers);
                else if (type == "RNA")
                    n2 = this.addHELMRNAs(plugin, ss, chain, origin, renamedmonomers);
                else if (type == "CHEM")
                    n2 = this.addChem(plugin, s, chain, origin, renamedmonomers);
                else if (type == "BLOB")
                    n2 = this.addBlob(plugin, s, chain, origin, renamedmonomers, e.tag);

                if (n2 > 0) {
                    n += n2;
                    origin.y += 4 * plugin.jsd.bondlength;
                }

                if (!scil.Utils.isNullOrEmpty(e.tag) && groups[sid] == null) {
                    ++gi;
                    var g = "G" + gi;
                    sections[2] += "|" + g + "(" + sid + ")";
                    groups[sid] = g;
                    groupannotations[g] = e.tag;
                }
            }
        }

        // hydrogenpairs
        var hydrogenpairs = [];
        var connections = [];
        var connatoms = {};

        s = sections[1];
        if (!scil.Utils.isNullOrEmpty(s)) {
            var ss = s == "" ? [] : this.split(s, '|');
            // RNA1,RNA1,1:R1-21:R2
            for (var i = 0; i < ss.length; ++i) {
                if (ss[i].indexOf("pair") > 0) {
                    hydrogenpairs.push(ss[i]);
                    continue;
                }

                var e = this.detachAnnotation(ss[i]);
                var c = this.parseConnection(e.str);
                if (c == null)
                    continue;

                if (isNaN(c.a1) && !/^G[0-9]+$/.test(c.chain1) && !scil.Utils.startswith(c.chain1, "BLOB")) {
                    // create group
                    if (groups[c.chain1] == null) {
                        ++gi;
                        var g = "G" + gi;
                        sections[2] += "|" + g + "(" + c.chain1 + ")";
                        groups[c.chain1] = g;
                    }
                    c.chain1 = groups[c.chain1];
                }
                if (isNaN(c.a2) && !/^G[0-9]+$/.test(c.chain2) && !scil.Utils.startswith(c.chain2, "BLOB")) {
                    // create group
                    if (groups[c.chain2] == null) {
                        ++gi;
                        var g = "G" + gi;
                        sections[2] += "|" + g + "(" + c.chain2 + ")";
                        groups[c.chain2] = g;
                    }
                    c.chain2 = groups[c.chain2];
                }

                c.tag = e.tag;
                connections.push(c);
                connatoms[c.chain1] = true;
                connatoms[c.chain2] = true;
            }
        }

        // groups, pairs, hydrogen bonds
        // RNA1,RNA2,2:pair-9:pair|RNA1,RNA2,5:pair-6:pair|RNA1,RNA2,8:pair-3:pair
        s = sections[2];
        var bondratios = [];
        if (!scil.Utils.isNullOrEmpty(s)) {
            var ss = s == "" ? [] : this.split(s, '|');
            for (var i = 0; i < ss.length; ++i) {
                if (scil.Utils.endswith(ss[i], ")") && /^[G|g][0-9]+[\(]/.test(ss[i])) {
                    // group
                    var p = ss[i].indexOf('(');
                    var c = ss[i].substr(p + 1, ss[i].length - p - 2);
                    var id = ss[i].substr(0, p);
                    if (connatoms[id] != null) {
                        var chain = this.createGroupForChains(plugin, chains, id, c);
                        if (chain != null)
                            chains[id] = chain;
                    }
                    else if (!this.parseBondRatios(bondratios, c)) { // bond ratio
                        // then group
                        var chain = this.createGroupForChains(plugin, chains, id, c, groupannotations[id]);
                        if (chain != null && chain.atoms.length == 1)
                            plugin.groupExpand(chain.atoms[0]);
                    }
                }
                else {
                    // pair
                    hydrogenpairs.push(ss[i]);
                }
            }
        }

        for (var i = 0; i < hydrogenpairs.length; ++i) {
            // pair
            var c = this.parseConnection(hydrogenpairs[i]);
            if (c == null || chains[c.chain1] == null || chains[c.chain2] == null || !scil.Utils.startswith(c.chain1, "RNA") || !scil.Utils.startswith(c.chain2, "RNA"))
                continue; //error

            var atom1 = chains[c.chain1].getAtomByAAID(c.a1);
            var atom2 = chains[c.chain2].getAtomByAAID(c.a2);
            if (atom1 == null || atom2 == null)
                continue; //error

            if (c.r1 != "pair" || c.r2 != "pair")
                continue; //error

            //chain.bonds.push(plugin.addBond(atom1, atom2, r1, r2));
            plugin.addHydrogenBond(atom1, atom2);
        }

        // connection
        for (var i = 0; i < connections.length; ++i) {
            var c = connections[i];
            if (c == null || chains[c.chain1] == null || chains[c.chain2] == null)
                continue; //error

            if (groups[c.chain1] != null)
                c.chain1 = groups[c.chain1];
            if (groups[c.chain2] != null)
                c.chain2 = groups[c.chain2];

            var chain1 = chains[c.chain1];
            var chain2 = chains[c.chain2];
            var atom1, atom2;
            var a1 = parseInt(c.a1);
            var a2 = parseInt(c.a2);
            if (a1 > 0 && !scil.Utils.startswith(c.chain1, "G") && !scil.Utils.startswith(c.chain1, "BLOB")) {
                atom1 = chain1.getAtomByAAID(c.a1);
            }
            else {
                atom1 = chain1.atoms[0];
                c.r1 = c.a1 + ":" + c.r1;
            }
            if (a2 > 0 && !scil.Utils.startswith(c.chain2, "G") && !scil.Utils.startswith(c.chain2, "BLOB")) {
                atom2 = chain2.getAtomByAAID(c.a2);
            }
            else {
                atom2 = chain2.atoms[0];
                c.r2 = c.a2 + ":" + c.r2;
            }
            if (atom1 == null || atom2 == null)
                continue; //error

            if (c.r1 == null || c.r2 == null)
                continue; // error
            var r1 = scil.Utils.startswith(c.r1, "R") ? parseInt(c.r1.substr(1)) : c.r1;
            var r2 = scil.Utils.startswith(c.r2, "R") ? parseInt(c.r2.substr(1)) : c.r2;
            var b = plugin.addBond(atom1, atom2, r1, r2);
            if (b != null) {
                b.tag = c.tag;
                var bondratio = this.findBondRatio(bondratios, groups, c.chain1, c.chain2);
                if (bondratio != null) {
                    b.ratio1 = bondratio.ratio1 != null ? bondratio.ratio1 : org.helm.webeditor.defaultbondratio;
                    b.ratio2 = bondratio.ratio2 != null ? bondratio.ratio2 : org.helm.webeditor.defaultbondratio;
                }
            }
        }

        // annotation
        s = sections[3];
        if (!scil.Utils.isNullOrEmpty(s)) {
            var ann = scil.Utils.eval(s);
            if (ann != null) {
                // HELM 2.0
                for (var k in ann) {
                    var chain = chains[k];
                    if (chain != null && chain.type == "RNA") {
                        var strandtype = ann[k].strandtype;
                        if (strandtype == "ss" || strandtype == "as")
                            chain.atoms[0].bio.annotation = "5'" + strandtype;
                    }
                }
            }
            else {
                // HELM 1.0
                var ss = this.split(s, '|');
                for (var i = 0; i < ss.length; ++i) {
                    var s = ss[i];
                    p = s.indexOf("{");
                    var chn = s.substr(0, p);
                    s = s.substr(p);
                    if (s == "{ss}" || s == "{as}") {
                        var chain = chains[chn];
                        if (chain != null && chain.type == "RNA")
                            chain.atoms[0].bio.annotation = "5'" + s.substr(1, s.length - 2);
                    }
                }
            }
        }

        return n;
    },

    findBondRatio: function (bondratios, groups, c1, c2) {
        for (var i = 0; i < bondratios.length; ++i) {
            var r = bondratios[i];
            var a1 = r.c1;
            var a2 = r.c2;
            if (groups[a1] != null)
                a1 = groups[a1];
            if (groups[a2] != null)
                a2 = groups[a2];

            if (a1 == c1 && a2 == c2)
                return { ratio1: r.ratio1, ratio2: r.ratio2 };
            else if (a1 == c2 && a2 == c1)
                return { ratio1: r.ratio2, ratio2: r.ratio1 };
        }
        return null;
    },

    parseBondRatios: function (bondratios, s) {
        var p = s.indexOf('+');
        if (p < 0)
            p = s.indexOf(',');
        if (p <= 0)
            return false;

        var ret = {};
        var s1 = s.substr(0, p);
        var s2 = s.substr(p + 1);

        p = s1.indexOf(':');
        if (p > 0) {
            ret.c1 = s1.substr(0, p);
            ret.ratio1 = s1.substr(p + 1);
        }
        else {
            ret.c1 = s1;
        }

        p = s2.indexOf(':');
        if (p > 0) {
            ret.c2 = s2.substr(0, p);
            ret.ratio2 = s2.substr(p + 1);
        }
        else {
            ret.c2 = s2;
        }

        bondratios.push(ret);
        return true;
    },

    createGroupForChains: function (plugin, chains, chainid, c, tag) {
        var logic = null;
        var ss = this.splitString(c, "+");
        if (ss.length > 1) {
            logic = "and";
        }
        else {
            ss = this.splitString(c, ",");
            if (ss.length > 1)
                logic = "or";
        }

        var allatoms = [];
        var atom = null;
        for (var i = 0; i < ss.length; ++i) {
            var ratio = null;

            var s = ss[i];
            var p = s.indexOf(':');
            if (p > 0) {
                ratio = s.substr(p + 1);
                s = s.substr(0, p);
            }

            var chain = chains[s];
            if (chain == null)
                continue; // error

            var atoms = [];
            for (var k = 0; k < chain.atoms.length; ++k) {
                atoms.push(chain.atoms[k]);
                if (chain.bases[k] != null)
                    atoms.push(chain.bases[k]);
            }

            //allatoms = allatoms.concat(atoms);
            var g2 = plugin.createGroup2(atoms, false);
            if (g2 != null) {
                g2.ratio = ratio;
                g2.tag = tag;

                var a2 = plugin.collapseGroup(g2);
                a2._aaid = 1;
                allatoms.push(a2);

                atom = a2;
            }
        }

        if (allatoms.length > 1) {
            var g = plugin.createGroup2(allatoms, false);
            if (g == null)
                return null;

            g.ratio = logic;
            g.tag = tag;
            atom = plugin.collapseGroup(g);
            atom._aaid = 1;
        }

        var chain = new org.helm.webeditor.Chain(ss[chainid]);
        chain.atoms.push(atom);
        return chain;
    },

    splitString: function (s, separators) {
        var ret = [];
        var w = "";
        for (var i = 0; i < s.length; ++i) {
            var c = s.substr(i, 1);
            if (separators.indexOf(c) >= 0) {
                ret.push(w);
                w = "";
            }
            else {
                w += c;
            }
        }

        if (ret.length == 0 || w.length > 0)
            ret.push(w);
        return ret;
    },

    /**
    * Split components (internal use)
    * @function split
    */
    split: function (s, sep) {
        var ret = [];
        // PEPTIDE1{G.[C[13C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|].T}$$$$

        var frag = "";
        var parentheses = 0;
        var bracket = 0;
        var braces = 0;
        var quote = 0;
        for (var i = 0; i < s.length; ++i) {
            var c = s.substr(i, 1);
            if (c == sep && bracket == 0 && parentheses == 0 && braces == 0 && quote == 0) {
                ret.push(frag);
                frag = "";
            }
            else {
                frag += c;
                if (quote > 0) {
                    if (c == '\\' && i + 1 < s.length) {
                        ++i;
                        var c2 = s.substr(i, 1);
                        frag += c2;
                        c += c2;
                    }
                }

                if (c == '\"') {
                    if (!(i > 0 && s.substr(i - 1, 1) == '\\'))
                        quote = quote == 0 ? 1 : 0;
                }
                else if (c == '[')
                    ++bracket;
                else if (c == ']')
                    --bracket;
                else if (c == '(')
                    ++parentheses;
                else if (c == ')')
                    --parentheses;
                else if (c == '{')
                    ++braces;
                else if (c == '}')
                    --braces;
            }
        }

        ret.push(frag);
        return ret;
    },

    /**
    * Parse HELM connection (internal use)
    * @function parseConnection
    */
    parseConnection: function (s) {
        var tt = s.split(',');
        if (tt.length != 3)
            return null; // error

        var tt2 = tt[2].split('-');
        if (tt2.length != 2)
            return null; // error

        var c1 = tt2[0].split(':');
        var c2 = tt2[1].split(':');
        if (c1.length != 2 || c2.length != 2)
            return null; // error

        return { chain1: tt[0], chain2: tt[1], a1: c1[0], r1: c1[1], a2: c2[0], r2: c2[1] };
        //return { chain1: tt[0], chain2: tt[1], a1: parseInt(c1[0]), r1: c1[1], a2: parseInt(c2[0]), r2: c2[1] };
    },

    /**
    * Split chars (internal use)
    * @function splitChars
    */
    splitChars: function (s, separator) {
        var ss = [];
        if (separator == null) {
            for (var i = 0; i < s.length; ++i)
                ss.push(s.substr(i, 1));
        }
        else {
            ss = s.split(separator);
        }
        return ss;
    },

    /**
    * Remove bracket (internal use)
    * @function trimBracket
    */
    trimBracket: function (s) {
        if (s != null && scil.Utils.startswith(s, "[") && scil.Utils.endswith(s, "]"))
            return s.substr(1, s.length - 2);
        return s;
    },

    /**
    * Make a renamed monomer (internal use)
    * @function getRenamedMonomer
    */
    getRenamedMonomer: function (type, elem, monomers) {
        if (monomers == null || monomers.length == 0)
            return elem;

        elem = org.helm.webeditor.IO.trimBracket(elem);
        for (var i = 0; i < monomers.length; ++i) {
            var m = monomers[i];
            if (m.oldname == elem)
                return m.id;
        }
        return elem;
    },

    /**
    * Remove annotation (internal use)
    * @function detachAnnotation
    */
    detachAnnotation: function (s) {
        var ret = this._detachAppendix(s, '\"');
        if (ret.tag != null)
            return ret;

        var r = this._detachAppendix(s, '\'');
        return { tag: ret.tag, repeat: r.tag, str: r.str };
    },

    _detachAppendix: function (s, c) {
        var tag = null;
        if (scil.Utils.endswith(s, c)) {
            var p = s.length - 1;
            while (p > 0) {
                p = s.lastIndexOf(c, p - 1);
                if (p <= 0 || s.substr(p - 1, 1) != '\\')
                    break;
            }

            if (p > 0 && p < s.length - 1) {
                tag = s.substr(p + 1, s.length - p - 2);
                s = s.substr(0, p);
            }
        }
        if (tag != null)
            tag = tag.replace(new RegExp("\\" + c, "g"), c);
        return { tag: this.unescape(tag), str: s };
    },

    unescape: function (s) {
        if (scil.Utils.isNullOrEmpty(s))
            return s;

        return s.replace(/[\\]./g, function (m) {
            switch (m) {
                case "\\r":
                    return "\r";
                case "\\n":
                    return "\n";
                case "\\t":
                    return "\t";
                default:
                    return m.substr(1);
            }
        });
    },

    escape: function (s) {
        if (scil.Utils.isNullOrEmpty(s))
            return s;

        return s.replace(/[\"|\'|\\|\r|\n|\t]/g, function (m) {
            switch (m) {
                case "\r":
                    return "\\r";
                case "\n":
                    return "\\n";
                case "\t":
                    return "\\t";
                default:
                    return "\\" + m;
            }
        });
    },

    /**
    * Add a monomer (internal use)
    * @function addNode
    */
    addNode: function (plugin, chain, atoms, p, type, elem, renamedmonomers) {
        var e = this.detachAnnotation(elem);
        a2 = plugin.addNode(p, type, this.getRenamedMonomer(type, e.str, renamedmonomers));
        if (a2 == null)
            throw "Failed to creating node: " + e.str;

        a2.tag = e.tag;
        atoms.push(a2);
        a2._aaid = chain.atoms.length + chain.bases.length;
        return a2;
    },

    /**
    * Add a CHEM node (internal use)
    * @function addChem
    */
    addChem: function (plugin, name, chain, origin, renamedmonomers) {
        this.addNode(plugin, chain, chain.atoms, origin.clone(), org.helm.webeditor.HELM.CHEM, name, renamedmonomers);
        return 1;
    },

    /**
    * Add a BLOB node (internal use)
    * @function addBlob
    */
    addBlob: function (plugin, name, chain, origin, renamedmonomers, annotation) {
        var e = this.detachAnnotation(name);
        var a = this.addNode(plugin, chain, chain.atoms, origin.clone(), org.helm.webeditor.HELM.BLOB, "Blob", renamedmonomers);
        a.bio.blobtype = e.str == "Blob" || e.str == "[Blob]" ? null : e.str;
        if (!scil.Utils.isNullOrEmpty(a.tag))
            a.tag = e.tag;
        else if (!scil.Utils.isNullOrEmpty(annotation))
            a.tag = annotation;
        return 1;
    },

    /**
    * Add Amino Acid (internal use)
    * @function addAAs
    */
    addAAs: function (plugin, ss, chain, origin, renamedmonomers) {
        var mol = plugin.jsd.m;
        var loop = { n: 0, firstatom: null, a1: null, a2: null, p: origin.clone(), delta: org.helm.webeditor.bondscale * plugin.jsd.bondlength };
        for (var i = 0; i < ss.length; ++i) {
            if (i == ss.length - 1 && ss[i] == ">") {
                if (loop.firstatom != loop.a1)
                    chain.bonds.push(plugin.addBond(loop.a1, loop.firstatom, 2, 1));
                break;
            }

            var e = this.detachAnnotation(ss[i]);

            var atoms = [];
            var rect = new JSDraw2.Rect(loop.p.x + loop.delta / 2, loop.p.y - loop.delta, 0, loop.delta * 2);
            if (scil.Utils.startswith(e.str, "(") && scil.Utils.endswith(e.str, ")")) {
                // dealing with repeat: PEPTIDE1{S.(D.F)'2-13'.A.S.D.F}$$$$V2.0

                // I#12364
                var ss2;
                var s = e.str.substr(1, e.str.length - 2);
                if (s.indexOf(',') > 0 || s.indexOf('+') > 0)
                    ss2 = [s]; // PEPTIDE1{A.A.A.A.(A:1.1,G:69.5,W:25.5,[Aha]:3.9)}$$$$
                else
                    ss2 = this.splitChars(s, '.'); // PEPTIDE1{A.A.A.A.(A.G)'2'.T}$$$$

                if (ss2.length == 1) {
                    atoms.push(this._addOneAA(plugin, chain, e.str, null, renamedmonomers, loop));
                }
                else {
                    for (var k = 0; k < ss2.length; ++k)
                        atoms.push(this._addOneAA(plugin, chain, ss2[k], null, renamedmonomers, loop));
                }
            }
            else {
                atoms.push(this._addOneAA(plugin, chain, e.str, e.tag, renamedmonomers, loop));
            }

            if (!scil.Utils.isNullOrEmpty(e.repeat)) {
                rect.width = loop.p.x + loop.delta / 2 - rect.left;
                var br = new JSDraw2.Bracket(null, rect);
                br.atoms = atoms;
                mol.addGraphics(br);
                br.createSubscript(mol, e.repeat);
            }
        }

        return loop.n;
    },

    _addOneAA: function (plugin, chain, s, tag, renamedmonomers, loop) {
        loop.p.x += loop.delta;
        var a = this.addNode(plugin, chain, chain.atoms, loop.p.clone(), org.helm.webeditor.HELM.AA, s, renamedmonomers);
        loop.a2 = a;
        loop.a2.tag = tag;

        if (loop.a1 != null)
            chain.bonds.push(plugin.addBond(loop.a1, loop.a2, 2, 1));

        if (loop.firstatom == null)
            loop.firstatom = loop.a2;

        loop.a1 = loop.a2;
        loop.a1.bio.id = ++loop.n;
        return a;
    },

    /**
    * Add RNA HELM string (internal use)
    * @function addHELMRNAs
    */
    addHELMRNAs: function (plugin, ss, chain, origin, renamedmonomers) {
        var mol = plugin.jsd.m;
        var loop = { n: 0, count: 0, firstatom: null, a1: null, a2: null, a3: null, p: origin.clone(), delta: org.helm.webeditor.bondscale * plugin.jsd.bondlength };
        for (var i = 0; i < ss.length; ++i) {
            var e = this.detachAnnotation(ss[i]);
            if (scil.Utils.startswith(e.str, "(") && scil.Utils.endswith(e.str, ")")) {
                var atoms = [];
                var rect = new JSDraw2.Rect(loop.p.x + loop.delta / 2, loop.p.y - loop.delta, 0, loop.delta * 3);

                var ss2 = this.splitChars(e.str.substr(1, e.str.length - 2), '.');
                for (var k = 0; k < ss2.length; ++k)
                    this._addOneHELMRNA(plugin, chain, ss2[k], renamedmonomers, loop, atoms);

                if (!scil.Utils.isNullOrEmpty(e.repeat)) {
                    rect.width = loop.p.x + loop.delta / 2 - rect.left;
                    var br = new JSDraw2.Bracket(null, rect);
                    br.atoms = atoms;
                    mol.addGraphics(br);
                    br.createSubscript(mol, e.repeat);
                }
            }
            else {
                this._addOneHELMRNA(plugin, chain, ss[i], renamedmonomers, loop, []);
            }
        }

        return loop.count;
    },

    _addOneHELMRNA: function (plugin, chain, s, renamedmonomers, loop, atoms) {
        var combo = this.splitCombo(s);
        for (var k = 0; k < combo.length; ++k) {
            var c = combo[k];
            var m = org.helm.webeditor.Monomers.getMonomer(org.helm.webeditor.HELM.SUGAR, c.symbol);
            if (m != null) {
                // sugar
                loop.p.x += loop.delta;
                loop.a2 = this.addNode(plugin, chain, chain.atoms, loop.p.clone(), org.helm.webeditor.HELM.SUGAR, c.symbol, renamedmonomers);
                if (loop.a1 != null)
                    chain.bonds.push(plugin.addBond(loop.a1, loop.a2, 2, 1));
                loop.a1 = loop.a2;

                if (!scil.Utils.isNullOrEmpty(c.base)) {
                    // base
                    loop.a3 = this.addNode(plugin, chain, chain.bases, org.helm.webeditor.Interface.createPoint(loop.p.x, loop.p.y + loop.delta), org.helm.webeditor.HELM.BASE, c.base, renamedmonomers);
                    plugin.addBond(loop.a1, loop.a3, 3, 1);
                    loop.a3.bio.id = ++loop.n;
                    ++loop.count;

                    atoms.push(loop.a3);
                }
            }
            else {
                if (!scil.Utils.isNullOrEmpty(c.base))
                    throw "Base attached to Linker: " + s;

                // linker
                var biotype = s == "*" ? org.helm.webeditor.HELM.NUCLEOTIDE : org.helm.webeditor.HELM.LINKER;
                loop.p.x += loop.delta;
                loop.a2 = this.addNode(plugin, chain, chain.atoms, loop.p.clone(), biotype, c.symbol, renamedmonomers);
                chain.bonds.push(plugin.addBond(loop.a1, loop.a2, 2, 1));
                loop.a1 = loop.a2;
                ++loop.count;
            }

            atoms.push(loop.a1);
            loop.a1.tag = c.tag;
        }
    },

    /**
    * Add RNA sequence (internal use)
    * @function addRNAs
    */
    addRNAs: function (plugin, ss, chain, origin, sugar, linker) {
        var n = 0;

        if (scil.Utils.isNullOrEmpty(sugar))
            sugar = "R";
        if (scil.Utils.isNullOrEmpty(linker) || linker == "null")
            linker = "P";

        var firstatom = null;
        var a1 = null;
        var a2 = null;
        var m = plugin.jsd.m;
        var delta = org.helm.webeditor.bondscale * plugin.jsd.bondlength;
        var p = origin.clone();
        for (var i = 0; i < ss.length; ++i) {
            if (i == ss.length - 1 && ss[i] == ">") {
                if (firstatom != a1) {
                    // linker
                    p.x += delta;
                    var a0 = this.addNode(plugin, chain, chain.atoms, p.clone(), org.helm.webeditor.HELM.LINKER, linker);
                    chain.bonds.push(plugin.addBond(a1, a0, 2, 1));

                    chain.bonds.push(plugin.addBond(a0, firstatom, 2, 1));
                }
                break;
            }

            // 1. linker
            if (a1 != null) {
                p.x += delta;
                var a0 = this.addNode(plugin, chain, chain.atoms, p.clone(), org.helm.webeditor.HELM.LINKER, linker);
                chain.bonds.push(plugin.addBond(a1, a0, 2, 1));
                a1 = a0;
            }

            // 2. sugar
            p.x += delta;
            a2 = this.addNode(plugin, chain, chain.atoms, p.clone(), org.helm.webeditor.HELM.SUGAR, sugar);
            if (a1 != null)
                chain.bonds.push(plugin.addBond(a1, a2, 2, 1));
            a1 = a2;

            if (firstatom == null)
                firstatom = a1;

            // 3. base
            a3 = this.addNode(plugin, chain, chain.bases, org.helm.webeditor.Interface.createPoint(p.x, p.y + delta), org.helm.webeditor.HELM.BASE, ss[i]);
            plugin.addBond(a1, a3, 3, 1);

            a3.bio.id = ++n;
        }

        return n;
    },

    /**
    * Split a RNA Combo (internal use)
    * @function splitCombo
    */
    splitCombo: function (s) {
        var ret = [];

        var m = null;
        var i = 0;
        while (i < s.length) {
            var c = s.substr(i, 1);
            if (c == '(') {
                if (i == 0)
                    throw "Invalid combo: " + s;
                var p = s.indexOf(')', i + 1);
                if (p <= i)
                    throw "Invalid combo: " + s;

                if (ret[ret.length - 1].base == null)
                    ret[ret.length - 1].base = s.substr(i + 1, p - i - 1);
                else
                    ret.push({ symbol: s.substr(i, p - i + 1) });

                i = p;
            }
            else if (c == '[') {
                var p = s.indexOf(']', i + 1);
                if (p <= i)
                    throw "Invalid combo: " + s;
                ret.push({ symbol: s.substr(i, p - i + 1) });
                i = p;
            }
            else if (c == '\"') {
                var p = s.indexOf('\"', i + 1);
                if (p <= i)
                    throw "Invalid combo: " + s;
                ret[ret.length - 1].tag = s.substr(i + 1, p - i - 1);
                i = p;
            }
            else {
                ret.push({ symbol: c });
            }

            ++i;
        }

        return ret;
    },

    /**
    * Compress a string using Pako (internal use)
    * @function compressGz
    */
    compressGz: function (s) {
        if (scil.Utils.isNullOrEmpty(s))
            return null;

        if (typeof pako != "undefined") {
            try {
                var buf = pako.deflate(s, { gzip: true });
                return btoa(String.fromCharCode.apply(null, buf));
            }
            catch (e) {
            }
        }
        return null;
    },

    /**
    * Decompress a string using pako (internal use)
    * @function uncompressGz
    */
    uncompressGz: function (b64Data) {
        if (scil.Utils.isNullOrEmpty(b64Data))
            return null;

        if (typeof pako == "undefined")
            return null;

        try {
            var strData = atob(b64Data);
            var charData = strData.split('').map(function (x) { return x.charCodeAt(0); });
            var binData = new Uint8Array(charData);
            var data = pako.inflate(binData);
            return String.fromCharCode.apply(null, new Uint16Array(data));
        }
        catch (e) {
            return null;
        }
    }
};
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* MonomerExplorer class
* @class org.helm.webeditor.MonomerExplorer
*/
org.helm.webeditor.MonomerExplorer = scil.extend(scil._base, {
    /**
    * @constructor MonomerExplorer
    * @param {DOM} parent - The parent element
    * @param {Plugin} plugin - The Plugin object
    * @param {dict} options - options on how to render the Monomer Explorer
    */
    constructor: function (parent, plugin, options) {
        this.plugin = plugin;
        this.options = options == null ? {} : options;
        this.height = null;
        var w = this.options.monomerwidth > 0 ? this.options.monomerwidth : 50;
        this.kStyle = { borderRadius: "5px", border: "solid 1px gray", backgroundRepeat: "no-repeat", display: "table", width: w, height: w, float: "left", margin: 2 };

        if (this.options.mexuseshape)
            this.kStyle.border = null;

        //this.lastselect = {};
        this.selected = {};
        this.selected[org.helm.webeditor.HELM.BASE] = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.BASE);
        this.selected[org.helm.webeditor.HELM.LINKER] = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.LINKER);
        this.selected[org.helm.webeditor.HELM.SUGAR] = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.SUGAR);
        this.selected[org.helm.webeditor.HELM.AA] = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.AA);
        this.selected[org.helm.webeditor.HELM.CHEM] = org.helm.webeditor.Monomers.getDefaultMonomer(org.helm.webeditor.HELM.CHEM);

        var me = this;
        this.div = scil.Utils.createElement(parent, "div", null, { fontSize: this.options.mexfontsize == null ? "90%" : this.options.mexfontsize });
        if (this.options.mexfind) {
            var d = scil.Utils.createElement(this.div, "div", null, { background: "#eee", borderBottom: "solid 1px gray", padding: "4px 0 4px 0" });
            var tbody = scil.Utils.createTable(d, 0, 0);
            var tr = scil.Utils.createElement(tbody, "tr");
            scil.Utils.createElement(tr, "td", "Quick Replace:", null, { colSpan: 3 });
            this.findtype = scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "select", null, { width: 100 });
            scil.Utils.listOptions(this.findtype, org.helm.webeditor.monomerTypeList(), null, true, false);

            tr = scil.Utils.createElement(tbody, "tr");
            this.findinput = scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "input", null, { width: 60 });
            scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "span", "&rarr;");
            this.findreplace = scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "input", null, { width: 60 });
            scil.Utils.createButton(scil.Utils.createElement(tr, "td", null, { textAlign: "right" }), { label: "Update", onclick: function () { me.findReplace(); } });
        }
        if (this.options.mexfilter != false) {
            var d = scil.Utils.createElement(this.div, "div", null, { background: "#eee", borderBottom: "solid 1px gray", padding: "4px 0 4px 0" });
            var tbody = scil.Utils.createTable(d, 0, 0);
            var tr = scil.Utils.createElement(tbody, "tr");
            scil.Utils.createElement(tr, "td", JSDraw2.Language.res("Filter") + ":", { paddingLeft: "5px" });
            this.filterInput = scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "input");
            scil.connect(this.filterInput, "onkeyup", function (e) { me.filter(e); });
        }

        var tabs = [];
        if (this.options.mexmonomerstab)
            tabs.push({ caption: "Monomers", tabkey: "monomers" });
        else
            this.addMonomerTabs(tabs);
        tabs.push({ caption: "Rules", tabkey: "rule" });

        var width = this.options.width != null ? this.options.width : 300;
        this.height = this.options.height != null ? this.options.height : 400;
        this.tabs = new scil.Tabs(scil.Utils.createElement(this.div, "div", null, { padding: "5px" }), {
            onShowTab: function (td) { me.onShowTab(td); },
            tabpadding: this.options.mexmonomerstab ? "10px" : "5px 2px 1px 2px",
            tabs: tabs,
            marginBottom: 0
        });

        this.dnd = this.createDnD(this.div);
        scil.connect(document.body, "onmousemove", function (e) { me.showMol(e); });

        org.helm.webeditor.MonomerExplorer.loadNucleotides();
    },

    /**
    * Add Monomer Tabs (internal use)
    * @function addMonomerTabs
    */
    addMonomerTabs: function (tabs) {
        if (this.options.mexfavoritetab != false)
            tabs.push({ caption: "Favorite", tabkey: "favorite" });

        tabs.push({ caption: "Chem", tabkey: "chem" });
        tabs.push({ caption: "Peptide", tabkey: "aa" });
        tabs.push({ caption: "RNA", tabkey: "rna" });
    },

    /**
    * Find and replace monomer (internal use)
    * @function findReplace
    */
    findReplace: function () {
        this.plugin.replaceAll(this.findinput.value, this.findreplace.value, this.findtype.value);
    },

    /**
    * Filter monomers (internal use)
    * @function filter
    */
    filter: function (e) {
        var key = this.tabs.currentTabKey();
        if (key == "rule") {
            org.helm.webeditor.RuleSet.filterRules(this.rules, this.filterInput.value, this.rules_category.value);
        }
        else {
            this.filterGroup(this.filterInput.value);
        }
    },

    /**
    * Filter a group (internal use)
    * @function filterGroup
    */
    filterGroup: function (s) {
        if (s == "")
            s = null;

        var groups = this.curtab.clientarea.className == "filtergroup" ? [this.curtab.clientarea] : this.curtab.clientarea.getElementsByClassName("filtergroup");
        for (var k = 0; k < groups.length; ++k) {
            var startingwith = [];
            var containing = [];
            var hidden = [];

            var parent = groups[k];
            for (var i = 0; i < parent.childNodes.length; ++i) {
                var d = parent.childNodes[i];
                var name = scil.Utils.getInnerText(d);
                var html = scil.Utils.isNullOrEmpty(name) ? d.innerHTML : null;
                var f = 1;
                if (s != null) {
                    f = 0;
                    if (scil.Utils.startswith(name.toLowerCase(), s)) {
                        f = 1;
                    }
                    else if (name.toLowerCase().indexOf(s) >= 0) {
                        f = 2;
                    }
                    else if (s.length >= 3 || org.helm.webeditor.MonomerExplorer.filtername) {
                        var type = d.getAttribute("helm");
                        var set = type == org.helm.webeditor.MonomerExplorer.kNucleotide ? org.helm.webeditor.MonomerExplorer.nucleotides : org.helm.webeditor.Monomers.getMonomerSet(type);
                        var m = set[scil.helm.symbolCase(name)];
                        if (m != null && m.n != null) {
                            if (scil.Utils.startswith(m.n.toLowerCase(), s))
                                f = 1;
                            else if (m.n.toLowerCase().indexOf(s) >= 0)
                                f = 2;
                        }
                    }
                }

                if (f == 1)
                    startingwith.push({ id: name, div: d, html: html });
                else if (f == 2)
                    containing.push({ id: name, div: d, html: html });
                else
                    hidden.push(d);
            }

            startingwith.sort(org.helm.webeditor.MonomerExplorer.compareMonomers);
            if (containing.length > 0) {
                containing.sort(org.helm.webeditor.MonomerExplorer.compareMonomers);
                startingwith = startingwith.concat(containing);
            }

            var last = null;
            for (var i = 0; i < startingwith.length; ++i) {
                var d = startingwith[i];
                parent.insertBefore(d.div, parent.childNodes[i]);
                last = d.div;
                if (s != null)
                    d.div.firstChild.firstChild.innerHTML = this.highlightString(d.id, s);
                else
                    d.div.firstChild.firstChild.innerHTML = d.html != null ? d.html : d.id;
                d.div.style.display = "table";
            }

            for (var i = 0; i < hidden.length; ++i)
                hidden[i].style.display = "none";
        }
    },

    /**
    * Highlight search string (internal use)
    * @function highlightString
    */
    highlightString: function (s, q) {
        var p = s.toLowerCase().indexOf(q);
        if (p < 0)
            return s;

        return s.substr(0, p) + "<span style='background:yellow'>" + s.substr(p, q.length) + "</span>" + s.substr(p + q.length);
    },

    /**
    * Reload a tab (internal use)
    * @function reloadTab
    */
    reloadTab: function (type) {
        var key = null;
        switch (type) {
            case "nucleotide":
                key = type;
                break;
            case org.helm.webeditor.HELM.AA:
                key = "aa";
                break;
            case org.helm.webeditor.HELM.CHEM:
                key = "chem";
                break;
            case org.helm.webeditor.HELM.BASE:
                key = "base";
                break;
            case org.helm.webeditor.HELM.LINKER:
                key = "linker";
                break;
            case org.helm.webeditor.HELM.SUGAR:
                key = "sugar";
                break;
            default:
                return;
        }

        var td = this.tabs.findTab(key);
        if (td == null && this.monomerstabs != null)
            td = this.monomerstabs.findTab(key);
        if (td == null)
            td = this.rnatabs.findTab(key);

        if (td != null)
            this.onShowTab(td, true);
    },

    /**
    * Reload all tabs (internal use)
    * @function reloadTabs
    */
    reloadTabs: function () {
        var list = this.tabs.tr.childNodes;
        for (var i = 0; i < list.length; ++i) {
            var td = list[i];
            scil.Utils.removeAll(td.clientarea);
            td._childrencreated = false;
        }

        this.onShowTab(this.tabs.currenttab);
    },

    /**
    * resize Monomer Explorer (internal use)
    * @function resize
    */
    resize: function (height) {
        this.height = height;

        if (this.divRule != null)
            this.divRule.style.height = this.getHeight("rule") + "px";
        if (this.divFavorite != null)
            this.divFavorite.style.height = this.getHeight("favorite") + "px";
        if (this.divChem != null)
            this.divChem.style.height = this.getHeight("chem") + "px";
        if (this.divAA != null)
            this.divAA.style.height = this.getHeight("aa") + "px";

        if (this.rnatabs != null)
            this.rnatabs.resizeClientarea(0, this.getHeight("RNA"));
    },

    /**
    * Get the height of the Monomer Explorer (internal use)
    * @function getHeight
    */
    getHeight: function (key) {
        var d1 = this.options.mexmonomerstab ? 0 : 14;
        var d2 = this.options.mexmonomerstab ? 0 : 47;
        var d3 = this.options.mexmonomerstab ? 0 : 46;
        switch (key) {
            case "rule":
                return this.height - 19 + d1;
            case "favorite":
                return this.height - 33 + d2;
            case "chem":
                return this.height - 33 + d2;
            case "aa":
                return this.height - 33 + d2;
            case "RNA":
                return this.height - 59 + d3;
        }

        return this.height;
    },

    /**
    * Event handler when showing a tab (internal use)
    * @function onShowTab
    */
    onShowTab: function (td, forcerecreate) {
        if (td == null)
            return;

        this.filterInput.value = "";
        this.curtab = td;
        this.filterGroup("");

        var key = td.getAttribute("key");
        if (forcerecreate || key == "favorite" && org.helm.webeditor.MonomerExplorer.favorites.changed) {
            td._childrencreated = false;
            if (key == "favorite")
                org.helm.webeditor.MonomerExplorer.favorites.changed = false;
        }

        if (this.plugin != null && this.plugin.jsd != null)
            this.plugin.jsd.doCmd("helm_" + key);
        if (td._childrencreated)
            return;
        td._childrencreated = true;

        var me = this;
        var div = td.clientarea;
        scil.Utils.unselectable(div);
        scil.Utils.removeAll(div);

        if (key == "favorite") {
            this.divFavorite = scil.Utils.createElement(div, "div", null, { width: "100%", height: this.getHeight(key), overflowY: "scroll" });
            this.recreateFavorites(this.divFavorite);
        }
        else if (key == "rna") {
            var d = scil.Utils.createElement(div, "div");
            this.createMonomerGroup3(d, "RNA", 0, false);
        }
        else if (key == "nucleotide") {
            var dict = org.helm.webeditor.MonomerExplorer.loadNucleotides();
            var list = scil.Utils.getDictKeys(dict);
            this.createMonomerGroup4(div, key, list);
        }
        else if (key == "aa") {
            this.divAA = scil.Utils.createElement(div, "div", null, { width: "100%", height: this.getHeight(key), overflowY: "scroll" });
            dojo.connect(this.divAA, "onmousedown", function (e) { me.select(e); });
            dojo.connect(this.divAA, "ondblclick", function (e) { me.dblclick(e); });
            this.createMonomerGroup4(this.divAA, org.helm.webeditor.HELM.AA, null, false, this.options.mexgroupanalogs != false);
        }
        else if (key == "chem") {
            this.divChem = scil.Utils.createElement(div, "div", null, { width: "100%", height: this.getHeight(key), overflowY: "scroll" });
            this.createMonomerGroup(this.divChem, org.helm.webeditor.HELM.CHEM);
        }
        else if (key == "base") {
            this.createMonomerGroup4(div, org.helm.webeditor.HELM.BASE, null, null, this.options.mexgroupanalogs != false);
        }
        else if (key == "sugar") {
            this.createMonomerGroup4(div, org.helm.webeditor.HELM.SUGAR, null);
        }
        else if (key == "linker") {
            this.createMonomerGroup4(div, org.helm.webeditor.HELM.LINKER, null, true);
        }
        else if (key == "rule") {
            var toolbar = scil.Utils.createElement(div, "div", null, { background: "#ccc" });
            scil.Utils.createElement(toolbar, "span", "Category:");
            this.rules_category = scil.Utils.createElement(toolbar, "select");
            scil.Utils.listOptions(this.rules_category, org.helm.webeditor.RuleManager.categories);
            var me = this;
            scil.connect(this.rules_category, "onchange", function () { org.helm.webeditor.RuleSet.filterRules(me.rules, me.filterInput.value, me.rules_category.value) });

            this.divRule = scil.Utils.createElement(div, "div", null, { width: "100%", height: this.getHeight(key), overflowY: "scroll" });
            this.listRules();
        }
        else if (key == "monomers") {
            var d = scil.Utils.createElement(div, "div", null, { paddingTop: "5px" });

            if (this.options.canvastoolbar == false) {
                var b = scil.Utils.createElement(d, "div", "<img src='" + scil.Utils.imgSrc("helm/arrow.png") + "' style='vertical-align:middle'>Mouse Pointer", { cursor: "pointer", padding: "2px", border: "solid 1px gray", margin: "5px" });
                scil.connect(b, "onclick", function () { me.plugin.jsd.doCmd("lasso"); });
            }

            var tabs = [];
            this.addMonomerTabs(tabs);
            this.monomerstabs = new scil.Tabs(d, {
                onShowTab: function (td) { me.onShowTab(td); },
                tabpadding: "5px 2px 1px 2px",
                tabs: tabs,
                marginBottom: 0
            });
        }
    },

    listRules: function () {
        var me = this;
        this.rules = org.helm.webeditor.RuleSet.listRules(this, function (script) { me.plugin.applyRule(script); }, function (scripts) { me.plugin.applyRules(scripts); });
    },

    /**
    * Get monomers by natural analog (internal use)
    * @function getMonomerDictGroupByAnalog
    */
    getMonomerDictGroupByAnalog: function (type) {
        var set = org.helm.webeditor.Monomers.getMonomerSet(type);
        //for (var k in set)
        //    set[k].id = k;

        var ret = {};
        var aa = type == org.helm.webeditor.HELM.AA;
        if (aa) {
            ret["C-Term"] = [];
            ret["N-Term"] = [];
        }

        for (var k in set) {
            var m = set[k];
            var na = m.na;
            if (aa) {
                if (m.at.R1 == null)
                    na = "N-Term";
                else if (m.at.R2 == null)
                    na = "C-Term";
            }
            if (scil.Utils.isNullOrEmpty(na))
                na = "X";
            if (ret[na] == null)
                ret[na] = [];
            ret[na].push(m);
        }

        for (var k in ret)
            ret[k] = this.getMonomerNames(ret[k]);

        return ret;
    },

    /**
    * Get monomer list of a monomer type (internal use)
    * @function getMonomerList
    */
    getMonomerList: function (list, type, addnull) {
        if (list != null) {
            list.sort();
            return list;
        }

        var set = org.helm.webeditor.Monomers.getMonomerSet(type);
        //for (var k in set)
        //    set[k].id = k;
        list = scil.Utils.getDictValues(set);
        return this.getMonomerNames(list, addnull);
    },

    /**
    * Get monomer names (internal use)
    * @function getMonomerNames
    */
    getMonomerNames: function (list, addnull) {
        var ret = [];
        //if (addnull)
        //    ret.push("null");

        list.sort(org.helm.webeditor.MonomerExplorer.compareMonomers);
        for (var i = 0; i < list.length; ++i)
            ret.push(list[i].id);

        return ret;
    },

    /**
    * Create a monomer group (internal use)
    * @function createMonomerGroup
    */
    createMonomerGroup: function (div, type, list, addnull) {
        var me = this;
        list = this.getMonomerList(list, type, addnull);

        if (org.helm.webeditor.ambiguity) {
            if (type == org.helm.webeditor.HELM.CHEM)
                list.splice(0, 0, '*');
        }

        div.style.overflowY = "scroll";
        this._listMonomers(div, list, type, this.options.mexfavoritefirst);
        dojo.connect(div, "onmousedown", function (e) { me.select(e); });
        dojo.connect(div, "ondblclick", function (e) { me.dblclick(e); });
    },

    /**
    * inner loop creating a monomer group (internal use)
    * @function createMonomerGroup3
    */
    createMonomerGroup3: function (div, group, i, createbar) {
        var me = this;
        var parent = scil.Utils.createElement(div, "div");
        if (createbar) {
            var bar = scil.Utils.createElement(parent, "div", group + ":", { background: "#ddd", borderTop: "solid 1px #aaa", marginTop: i == 0 ? null : "1px" });
            if (i > 0)
                new scil.Resizable(bar, { direction: "y", mouseovercolor: "#aaf", onresize: function (delta, resizable) { return me.onresize(delta, i); } });
        }

        var d = scil.Utils.createElement(parent, "div");
        dojo.connect(d, "onmousedown", function (e) { me.select(e); });
        dojo.connect(d, "ondblclick", function (e) { me.dblclick(e); });

        if (group == "RNA") {
            var base = org.helm.webeditor.Monomers.bases["A"] == null ? "a" : "A";
            var linker = org.helm.webeditor.Monomers.linkers["P"] == null ? "p" : "P";
            var sugar = org.helm.webeditor.Monomers.sugars["R"] == null ? "r" : "R";

            var tabs = [
                    { caption: this.createRNATabCaption("nucleotide", "R(A)P"), tabkey: "nucleotide", onmenu: this.options.mexrnapinontab ? function (e) { me.onPinMenu(e); } : null },
                    { caption: this.createRNATabCaption("base", base), tabkey: "base" },
                    { caption: this.createRNATabCaption("sugar", sugar), tabkey: "sugar" },
                    { caption: this.createRNATabCaption("linker", linker), tabkey: "linker" }
                ];
            this.rnatabs = new scil.Tabs(scil.Utils.createElement(d, "div", null, { paddingTop: "5px" }), {
                onShowTab: function (td) { me.onShowTab(td); }, //function (td) { me.onShowRNATab(td); },
                tabpadding: "2px",
                tabs: tabs,
                marginBottom: 0,
                clientareaheight: this.getHeight("RNA")
            });
        }
        //else if (group == "Chem") {
        //    d.style.overflowY = "scroll";
        //    d.style.height = height + "px";
        //    var list = this.getMonomerList(null, org.helm.webeditor.HELM.CHEM);
        //    this._listMonomers(d, list, org.helm.webeditor.HELM.CHEM, true);
        //}
        //else if (group == "Peptide") {
        //    d.style.overflowY = "scroll";
        //    d.style.height = height + "px";
        //    this.createMonomerGroup4(d, org.helm.webeditor.HELM.AA, null, false, this.options.mexgroupanalogs != false);
        //    //var list = this.getMonomerList(null, org.helm.webeditor.HELM.AA);
        //    //this._listMonomers(d, list, org.helm.webeditor.HELM.AA, true);
        //}
    },

    /**
    * Create RNA Tab caption (internal use)
    * @function 
    */
    createRNATabCaption: function (type, label) {
        var half = " style='font-size: 80%;padding-left:20px;background-repeat:no-repeat;background-position:left center;background-image:";
        return "<div title='Nucleotide (Combined)' " + half + scil.Utils.imgSrc("img/helm_" + type.toLowerCase() + ".gif", true) + "'>" + label + "</div>";
    },

    /**
    * Event handler when pinning a Combo (internal use)
    * @function onPinMenu
    */
    onPinMenu: function (e) {
        if (this.pinmenu == null) {
            var me = this;
            var items = [{ caption: "Pin This Nucleotide"}];
            this.pinmenu = new scil.ContextMenu(items, function () { me.addNucleotide(); });
        }
        this.pinmenu.show(e.clientX, e.clientY);
    },

    /**
    * Inner loop creating a monomer group (internal use)
    * @function createMonomerGroup4
    */
    createMonomerGroup4: function (div, type, list, addnull, groupbyanalog) {
        if (groupbyanalog) {
            var dict = this.getMonomerDictGroupByAnalog(type);

            if (org.helm.webeditor.ambiguity) {
                if (type == org.helm.webeditor.HELM.AA)
                    dict['?'] = ['*', '_', 'X'];
                else if (type == org.helm.webeditor.HELM.BASE)
                    dict['?'] = ['_', 'N', '*'];
            }

            var list = [];
            if (this.options.mexfavoritefirst) {
                for (var k in dict) {
                    var list2 = dict[k];
                    for (var i = 0; i < list2.length; ++i) {
                        var a = list2[i];
                        if (org.helm.webeditor.MonomerExplorer.favorites.contains(a, type))
                            list.push(a);
                    }
                }
                this._listMonomer2(div, scil.Utils.imgTag("star.png"), list, type, 20);
            }

            list = scil.Utils.getDictKeys(dict);
            list.sort();
            var list2 = [];
            for (var i = 0; i < list.length; ++i) {
                var k = list[i];
                if (k == "C-Term" || k == "N-Term") {
                    list2.push(k);
                    continue;
                }
                this._listMonomer2(div, k, dict[k], type, 20);
            }

            for (var i = 0; i < list2.length; ++i) {
                var k = list2[i];
                this._listMonomer2(div, k, dict[k], type, 60);
            }
        }
        else {
            if (type == "nucleotide" && !this.options.mexrnapinontab) {
                var me = this;
                var d = this.createMonomerDiv(div, scil.Utils.imgTag("pin.png"), null, null, false);
                d.setAttribute("title", "Pin This Nucleotide");
                scil.connect(d, "onclick", function () { me.addNucleotide(); })
            }
            var list = this.getMonomerList(list, type, addnull);
            if (org.helm.webeditor.ambiguity) {
                if (type == org.helm.webeditor.HELM.SUGAR)
                    list.splice(0, 0, '*');
                else if (type == org.helm.webeditor.HELM.LINKER)
                    list.splice(0, 0, '*');
                else if (type == "nucleotide")
                    list.splice(0, 0, '*');
            }

            this._listMonomers(div, list, type, this.options.mexfavoritefirst);
        }
    },

    /**
    * Add a nucleotide (internal use)
    * @function addNucleotide
    */
    addNucleotide: function (tab) {
        var notation = this.getCombo();
        var dict = org.helm.webeditor.MonomerExplorer.nucleotides;
        for (var k in dict) {
            if (notation == dict[k]) {
                scil.Utils.alert("There is a defined nucleotide called: " + k);
                return;
            }
        }

        var me = this;
        scil.Utils.prompt2({
            caption: "Pin Nucleotide: " + notation,
            message: "Please give a short name for the nucleotide, " + notation,
            callback: function (s) { if (org.helm.webeditor.MonomerExplorer.addCustomNucleotide(s, notation)) me.reloadTab("nucleotide"); }
        });
    },

    /**
    * Inner loop listing all monomer of a monomer type (internal use)
    * @function 
    */
    _listMonomer2: function (div, k, list, type, width) {
        if (list.length == 0)
            return;

        var tbody = scil.Utils.createTable(div, 0, 0);
        var tr = scil.Utils.createElement(tbody, "tr");
        var left = scil.Utils.createElement(tr, "td", null, { verticalAlign: "top" });
        var right = scil.Utils.createElement(tr, "td", null, { verticalAlign: "top" });
        scil.Utils.createElement(left, "div", k, { width: width, background: "#eee", border: "solid 1px #aaa", textAlign: "center" });
        this._listMonomers(right, list, type);
    },

    /**
    * Create favorite monomer group (internal use)
    * @function createMonomerGroupFav
    */
    createMonomerGroupFav: function (div, caption, type) {
        var list = org.helm.webeditor.MonomerExplorer.favorites.getList(type);
        if (list == null || list.length == 0)
            return;

        list.sort();
        scil.Utils.createElement(div, "div", caption + ":", { background: "#ddd", border: "solid 1px #ddd" });
        var d = scil.Utils.createElement(div, "div", null, { display: "table", paddingBottom: "10px" });
        this._listMonomers(d, list, type, false);

        var me = this;
        dojo.connect(d, "onmousedown", function (e) { me.select(e); });
        dojo.connect(d, "ondblclick", function (e) { me.dblclick(e); });
    },

    /**
    * List a monomer group (internal use)
    * @function _listMonomers
    */
    _listMonomers: function (div, list, type, mexfavoritefirst) {
        div.className = "filtergroup";

        if (mexfavoritefirst) {
            var list2 = [];
            for (var i = 0; i < list.length; ++i) {
                if (org.helm.webeditor.MonomerExplorer.favorites.contains(list[i], type))
                    this.createMonomerDiv(div, list[i], type);
                else
                    list2.push(list[i]);
            }

            for (var i = 0; i < list2.length; ++i)
                this.createMonomerDiv(div, list2[i], type);
        }
        else {
            for (var i = 0; i < list.length; ++i)
                this.createMonomerDiv(div, list[i], type);
        }
    },

    /**
    * Recreate favorite monomers (internal use)
    * @function 
    */
    recreateFavorites: function (d) {
        this.createMonomerGroupFav(d, "Nucleotide", org.helm.webeditor.MonomerExplorer.kNucleotide);
        this.createMonomerGroupFav(d, "Base", org.helm.webeditor.HELM.BASE);
        this.createMonomerGroupFav(d, "Sugar", org.helm.webeditor.HELM.SUGAR);
        this.createMonomerGroupFav(d, "Linker", org.helm.webeditor.HELM.LINKER);
        this.createMonomerGroupFav(d, "Chemistry", org.helm.webeditor.HELM.CHEM);
        this.createMonomerGroupFav(d, "Peptide", org.helm.webeditor.HELM.AA);
    },

    /**
    * Create a monomer block (internal use)
    * @function createMonomerDiv
    */
    createMonomerDiv: function (parent, name, type, style, star) {
        var fav = org.helm.webeditor.MonomerExplorer.favorites.contains(name, type);

        if (style == null)
            style = scil.clone(this.kStyle);
        else
            style = scil.apply(scil.clone(this.kStyle), style);

        if (this.options.mexusecolor != false) {
            var color;
            var custom = org.helm.webeditor.MonomerExplorer.customnucleotides;
            if (type == "nucleotide" && custom != null && custom[name] != null)
                color = { backgroundcolor: "#afa" };
            else
                color = style.backgroundColor = org.helm.webeditor.Monomers.getColor2(type, name);
            style.backgroundColor = color == null ? null : color.backgroundcolor;
        }

        if (star != false)
            style.backgroundImage = scil.Utils.imgSrc("img/star" + (fav ? "" : "0") + ".png", true);

        var div = scil.Utils.createElement(parent, "div", null, style, { helm: type, bkcolor: style.backgroundColor, star: (star ? 1 : null) });
        scil.Utils.unselectable(div);

        if (this.options.mexuseshape)
            this.setMonomerBackground(div, 0);

        var d = scil.Utils.createElement(div, "div", null, { display: "table-cell", textAlign: "center", verticalAlign: "middle" });
        scil.Utils.createElement(d, "div", name, { overflow: "hidden", width: this.kStyle.width });

        return div;
    },

    /**
    * Set monomer block background (internal use)
    * @function setMonomerBackground
    */
    setMonomerBackground: function (div, f) {
        var type = div.getAttribute("helm");
        if (scil.Utils.isNullOrEmpty(type))
            return;

        var bk = type.toLowerCase();
        if (type != org.helm.webeditor.MonomerExplorer.kNucleotide)
            bk = bk.substr(bk.indexOf('_') + 1);
        div.style.backgroundImage = scil.Utils.imgSrc("img/mon-" + bk + f + ".png", true);
    },

    /**
    * Get the monomer block (internal use)
    * @function getMonomerDiv
    */
    getMonomerDiv: function (e) {
        var div = e.target || e.srcElement;
        if (div == null || div.tagName == null)
            return;

        for (var i = 0; i < 3; ++i) {
            var type = div.getAttribute("helm");
            if (!scil.Utils.isNullOrEmpty(type))
                break;
            div = div.tagName == "BODY" ? null : div.parentNode;
            if (div == null)
                break;
        }
        return scil.Utils.isNullOrEmpty(type) ? null : div;
    },

    /**
    * Enable DnD (internal use)
    * @function createDnD
    */
    createDnD: function (div) {
        var me = this;
        return new scil.DnD(div, {
            onstartdrag: function (e, dnd) {
                return me.getMonomerDiv(e);
            },
            oncreatecopy: function (e, dnd) {
                if (me.dnd.floatingbox == null) {
                    var maxZindex = scil.Utils.getMaxZindex();
                    var style = {
                        float: null, backgroundImage: null,
                        filter: 'alpha(opacity=80)', opacity: 0.8, color: org.helm.webeditor.MonomerExplorer.color,
                        backgroundColor: org.helm.webeditor.MonomerExplorer.backgroundcolor,
                        zIndex: (maxZindex > 0 ? maxZindex : 100) + 1, position: "absolute"
                    };
                    if (me.options.useshape)
                        style.backgroundColor = null;
                    me.dnd.floatingbox = me.createMonomerDiv(document.body, null, null, style, false);
                }
                me.dnd.floatingbox.style.display = "table";
                me.dnd.floatingbox.style.backgroundColor = org.helm.webeditor.MonomerExplorer.backgroundcolor;
                me.dnd.floatingbox.innerHTML = dnd.src.innerHTML;
                me.dnd.floatingbox.setAttribute("helm", dnd.src.getAttribute("helm"));
                if (me.options.useshape)
                    me.setMonomerBackground(me.dnd.floatingbox, 1);
                return me.dnd.floatingbox;

            },
            ondrop: function (e, dnd) {
                if (me.dnd.floatingbox == null)
                    return;

                me.dnd.floatingbox.style.display = "none";

                var src = e.target || e.srcElement;
                if (!scil.Utils.hasAnsestor(src, me.plugin.jsd.div))
                    return;

                var type = me.dnd.floatingbox.getAttribute("helm");
                me.plugin.dropMonomer(type, scil.Utils.getInnerText(me.dnd.floatingbox), e);
            },
            oncancel: function (dnd) {
                if (me.dnd.floatingbox == null)
                    return;

                me.dnd.floatingbox.style.display = "none";
                var type = me.dnd.floatingbox.getAttribute("helm");
            }
        });
    },

    /**
    * Show structure popup (internal use)
    * @function showMol
    */
    showMol: function (e) {
        var src = this.getMonomerDiv(e);
        if (src != null && !this.dnd.isDragging()) {
            var type = src.getAttribute("helm");
            var set = type == org.helm.webeditor.MonomerExplorer.kNucleotide ? org.helm.webeditor.MonomerExplorer.nucleotides : org.helm.webeditor.Monomers.getMonomerSet(type);
            var s = scil.Utils.getInnerText(src);
            var m = set[scil.helm.symbolCase(s)];
            if (m == null)
                m = set[s];
            org.helm.webeditor.MolViewer.show(e, type, m, s);
        }
        else {
            var src = e.srcElement || e.target;
            if (!scil.Utils.isChildOf(src, this.plugin.jsd.div))
                org.helm.webeditor.MolViewer.hide();
        }
    },

    /**
    * Tool to split a list (internal use)
    * @function splitLists
    */
    splitLists: function (set) {
        var lists = [[], [], [], []];
        for (var k in set) {
            var m = set[k];
            if (m.at.R1 == null)
                lists[2].push(k);
            else if (m.at.R2 == null)
                lists[3].push(k);
            else if (k.length == 1)
                lists[0].push(k);
            else
                lists[1].push(k);
        }

        return lists;
    },

    /**
    * Change favorites (internal use)
    * @function changeFavorite
    */
    changeFavorite: function (div) {
        var f = div.getAttribute("star") != "1";

        if (f) {
            div.setAttribute("star", "1");
            div.style.backgroundImage = scil.Utils.imgSrc("img/star.png", true);
        }
        else {
            div.setAttribute("star", "");
            div.style.backgroundImage = scil.Utils.imgSrc("img/star0.png", true);
        }

        var type = div.getAttribute("helm");
        var s = scil.Utils.getInnerText(div);
        org.helm.webeditor.MonomerExplorer.favorites.add(s, f, type);

        //this.reloadTab(type);
    },

    /**
    * Select a monomer (internal use)
    * @function select
    */
    select: function (e) {
        var div = this.getMonomerDiv(e);
        if (div != null) {
            var d = scil.Utils.getOffset(div, true);
            var scroll = scil.Utils.getParent(div.parentNode, "div");
            var dx = e.clientX - d.x + scroll.scrollLeft;
            var dy = e.clientY - d.y + scroll.scrollTop;
            if (dx >= 0 && dx < 16 && dy >= 0 && dy < 16) {
                // favorite
                this.changeFavorite(div);
                e.preventDefault();
                return;
            }
        }

        var helm = div == null ? null : div.getAttribute("helm");
        if (scil.Utils.isNullOrEmpty(helm))
            return;

        this.plugin.jsd.activate(true);

        var name = scil.Utils.getInnerText(div);
        if (helm == org.helm.webeditor.MonomerExplorer.kNucleotide) {
            if (name != "*") {
                var s = org.helm.webeditor.MonomerExplorer.nucleotides[name];
                var p1 = s.indexOf('(');
                var p2 = s.indexOf(")");
                var sugar = org.helm.webeditor.IO.trimBracket(s.substr(0, p1));
                var base = org.helm.webeditor.IO.trimBracket(s.substr(p1 + 1, p2 - p1 - 1));
                var linker = org.helm.webeditor.IO.trimBracket(s.substr(p2 + 1));

                if (scil.Utils.isNullOrEmpty(linker))
                    linker = "null";

                this.selected[org.helm.webeditor.HELM.BASE] = base;
                this.selected[org.helm.webeditor.HELM.LINKER] = linker;
                this.selected[org.helm.webeditor.HELM.SUGAR] = sugar;

                if (this.rnatabs != null) {
                    var tabs = this.rnatabs;
                    tabs.updateTabLabel("nucleotide", this.createRNATabCaption("nucleotide", s));
                    tabs.updateTabLabel("sugar", this.createRNATabCaption("sugar", sugar));
                    tabs.updateTabLabel("linker", this.createRNATabCaption("linker", linker));
                    tabs.updateTabLabel("base", this.createRNATabCaption("base", base));
                }
            }
        }
        else {
            name = org.helm.webeditor.IO.trimBracket(name);
            if (this.rnatabs != null) {
                var tab = null;
                var tabs = this.rnatabs;
                switch (helm) {
                    case org.helm.webeditor.HELM.SUGAR:
                        tab = "sugar";
                        break;
                    case org.helm.webeditor.HELM.LINKER:
                        tab = "linker";
                        break;
                    case org.helm.webeditor.HELM.BASE:
                        tab = "base";
                        break;
                }

                if (tab != null)
                    tabs.updateTabLabel(tab, this.createRNATabCaption(tab, name));
            }

            this.selected[helm] = name;
            if (tabs != null)
                tabs.updateTabLabel("nucleotide", this.createRNATabCaption("nucleotide", this.getCombo()));
        }

        if (this.lastdiv != null) {
            this.lastdiv.style.color = "";
            if (this.options.mexuseshape) {
                this.setMonomerBackground(this.lastdiv, 0);
            }
            else {
                var s = this.lastdiv.getAttribute("bkcolor");
                this.lastdiv.style.backgroundColor = s == null ? "" : s;
            }
        }
        if (this.options.mexuseshape)
            this.setMonomerBackground(div, 1);
        else
            div.style.backgroundColor = org.helm.webeditor.MonomerExplorer.backgroundcolor;
        div.style.color = org.helm.webeditor.MonomerExplorer.color;
        this.lastdiv = div;

        if (this.plugin != null && this.plugin.jsd != null) {
            switch (helm) {
                case org.helm.webeditor.HELM.AA:
                    this.plugin.jsd.doCmd("helm_aa");
                    break;
                case org.helm.webeditor.HELM.CHEM:
                    this.plugin.jsd.doCmd("helm_chem");
                    break;
                case org.helm.webeditor.HELM.BASE:
                    this.plugin.jsd.doCmd(this.options.alwaysdrawnucleotide ? "helm_nucleotide" : "helm_base");
                    break;
                case org.helm.webeditor.HELM.LINKER:
                    this.plugin.jsd.doCmd(this.options.alwaysdrawnucleotide ? "helm_nucleotide" : "helm_linker");
                    break;
                case org.helm.webeditor.HELM.SUGAR:
                    this.plugin.jsd.doCmd(this.options.alwaysdrawnucleotide ? "helm_nucleotide" : "helm_sugar");
                    break;
                case org.helm.webeditor.MonomerExplorer.kNucleotide:
                    this.plugin.jsd.doCmd("helm_nucleotide");
                    break;
            }
        }
    },

    /**
    * Get the RNA combo (internal use)
    * @function getCombo
    */
    getCombo: function () {
        var sugar = this.selected[org.helm.webeditor.HELM.SUGAR];
        var linker = this.selected[org.helm.webeditor.HELM.LINKER];
        var base = this.selected[org.helm.webeditor.HELM.BASE];
        var s = org.helm.webeditor.IO.getCode(sugar);
        if (!org.helm.webeditor.Monomers.hasR(org.helm.webeditor.HELM.SUGAR, sugar, "R3"))
            s += "()";
        else
            s += "(" + org.helm.webeditor.IO.getCode(base) + ")";
        if (linker != "null")
            s += org.helm.webeditor.IO.getCode(linker);
        return s;
    },

    /**
    * Event of double click (internal use)
    * @function dblclick
    */
    dblclick: function (e) {
        var div = this.getMonomerDiv(e);
        var helm = div == null ? null : div.getAttribute("helm");
        if (org.helm.webeditor.isHelmNode(helm)) {
            if (this.plugin.dblclickMomonor(helm, scil.Utils.getInnerText(div)) == 0)
                scil.Utils.beep();
        }
    }
});


scil.apply(org.helm.webeditor.MonomerExplorer, {
    kUseShape: false,
    kNucleotide: "nucleotide",
    backgroundcolor: "blue",
    color: "white",
    customnucleotides: null,
    filtername: false,
    favorites: new scil.Favorite("monomers", function (name, f, type) { org.helm.webeditor.MonomerExplorer.onAddFavorite(name, f, type); }),

    nucleotides: {
        A: "R(A)P",
        C: "R(C)P",
        G: "R(G)P",
        T: "R(T)P",
        U: "R(U)P"
    },

    /**
    * Compare if two monomers are same (internal use)
    * @function compareMonomers
    */
    compareMonomers: function (a, b) {
        if (a.id == b.id)
            return 0;
        else if (a.id == null)
            return -1;
        else if (b.id == null)
            return 1;
        else if (a.id.length != b.id.length && (a.id.length == 1 || b.id.length == 1))
            return a.id.length > b.id.length ? 1 : -1;
        else
            return scil.helm.symbolCase(a.id) > scil.helm.symbolCase(b.id) ? 1 : -1;
    },

    /**
    * Event handler on adding favorite (internal use)
    * @function onAddFavorite
    */
    onAddFavorite: function (name, f, type) {
        if (!f && type == "nucleotide" && this.customnucleotides != null && this.customnucleotides[name] != null) {
            delete this.customnucleotides[name];
            this.saveNucleotides();
        }
    },

    /**
    * Create a custom RNA combo (internal use)
    * @function addCustomNucleotide
    */
    addCustomNucleotide: function (name, notation) {
        name = scil.Utils.trim(name);
        if (scil.Utils.isNullOrEmpty(name)) {
            scil.Utils.alert("The short name cannot be blank");
            return false;
        }

        if (this.nucleotides[name] != null) {
            scil.Utils.alert("The short name is used for: " + this.nucleotides[name]);
            return false;
        }

        if (this.customnucleotides == null)
            this.customnucleotides = {};

        this.nucleotides[name] = notation;
        this.customnucleotides[name] = notation;
        this.saveNucleotides();
        this.favorites.add(name, true, "nucleotide");

        return true;
    },

    /**
    * Save custom RNA Combos(internal use)
    * @function saveNucleotides
    */
    saveNucleotides: function () {
        var s = scil.Utils.json2str(this.customnucleotides);
        scil.Utils.createCookie("scil_helm_nucleotides", s);
    },

    /**
    * Read all saved custom RNA Combo (internal use)
    * @function loadNucleotides
    */
    loadNucleotides: function () {
        if (this._nucleotidesloaded)
            return this.nucleotides;

        if (this.nucleotides == null)
            this.nucleotides = [];

        this._nucleotidesloaded = true;
        var s = scil.Utils.readCookie("scil_helm_nucleotides");
        this.customnucleotides = scil.Utils.eval(s);
        if (this.customnucleotides != null && this.customnucleotides.length == null) {
            var list = {};
            for (var k in this.customnucleotides) {
                if (this.nucleotides[k] == null) {
                    list[k] = this.customnucleotides[k];
                    this.nucleotides[k] = this.customnucleotides[k];
                }
            }
            this.customnucleotides = list;
        }
        return this.nucleotides;
    },

    /**
    * Show Monomer Explorer as a dialog (internal use)
    * @function showDlg
    */
    showDlg: function (jsd) {
        this.createDlg(jsd);
        this.dlg.show2({ owner: jsd, modal: false });
        jsd.helm.monomerexplorer = this.mex;
    },

    /**
    * Create Monomer Explorer dialog (internal use)
    * @function createDlg
    */
    createDlg: function (jsd) {
        if (this.dlg != null)
            return;

        var div = scil.Utils.createElement(null, "div", null, { width: 500 });
        this.dlg = new scil.Dialog("Monomer Explorer", div);
        this.dlg.show2({ owner: jsd, modal: false });

        this.mex = new org.helm.webeditor.MonomerExplorer(div, jsd.helm, { height: 350 });
        this.dlg.moveCenter();
    }
});
﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* MolViewer class
* @class org.helm.webeditor.MolViewer
*/
org.helm.webeditor.MolViewer = {
    dlg: null,
    jsd: null,
    molscale: 1,
    delay: 800,

    /**
    * Show structure view popup
    * @function show
    */
    show: function (e, type, m, code, ed, text) {
        this.clearTimer();
        var me = this;
        this.tm = setTimeout(function () { me.show2({ x: e.clientX, y: e.clientY }, type, m, code, ed, text); }, this.delay);
    },

    /**
    * Clear delay timer
    * @function clearTimer
    */
    clearTimer: function () {
        if (this.tm != null) {
            clearTimeout(this.tm);
            this.tm = null;
        }
    },

    /**
    * Inner implementation of display structure dialog (internal use)
    * @function show2
    */
    show2: function (xy, type, m, code, ed, a) {
        this.tm = null;
        if (m == null && a == null)
            return;

        if (ed != null && ed.contextmenu != null && ed.contextmenu.isVisible())
            return;

        this.create();

        if (this.cur != (type + "." + code) || !this.dlg.isVisible()) {
            this.cur = type + "." + code;

            if (m != null && typeof (m) == "string") {
                var s = m;
                m = { n: m, m: this.assemblyMol(s) };
            }

            var name = "";
            if (m != null) {
                name = m.n;
                if (name == null)
                    name = a.elem;
            }
            else {
                if (a.bio != null && !scil.Utils.isNullOrEmpty(a.bio.ambiguity))
                    name = a.bio.ambiguity;
            }

            var blobtype = "";
            if (a != null && type == org.helm.webeditor.HELM.BLOB && !scil.Utils.isNullOrEmpty(a.bio.blobtype))
                blobtype = "{" + a.bio.blobtype + "}";

            var fields = this.dlg.form.fields;
            this.dlg.show2({ title: "<div style='font-size:80%'>" + name + blobtype + "</div>", modal: false, immediately: true });

            var molfile = org.helm.webeditor.monomers.getMolfile(m);
            if (scil.Utils.isNullOrEmpty(molfile)) {
                fields.jsd.style.display = "none";
            }
            else {
                fields.jsd.style.display = "";
                fields.jsd.jsd.setXml(molfile);
            }

            var s = "";
            if (m != null && m.at != null) {
                for (var k in m.at)
                    s += "<tr><td>" + k + "=</td><td>&nbsp;" + m.at[k] + "</td></tr>";
            }

            if (s == "") {
                fields.rs.style.display = "none";
            }
            else {
                fields.rs.style.display = "";
                fields.rs.innerHTML = "<table cellspacing=0 cellpadding=0 style='font-size:80%'>" + s + "</table>";
            }

            var s = "";
            if (a != null) {
                if (!scil.Utils.isNullOrEmpty(a.tag))
                    s += "<div>" + a.tag + "</div>";
            }
            if (s == "") {
                fields.notes.style.display = "none";
            }
            else {
                fields.notes.style.display = "";
                fields.notes.innerHTML = s;

                fields.notes.style.borderTop = fields.rs.style.display == "" ? "solid 1px gray" : "";
            }
        }

        var scroll = scil.Utils.scrollOffset();
        this.dlg.moveTo(xy.x + scroll.x + 10, xy.y + scroll.y + 10);
    },

    /**
    * Assembly molecule (internal use)
    * @function assemblyMol
    */
    assemblyMol: function (s) {
        var p1 = s.indexOf('(');
        var p2 = s.indexOf(")");
        var sugar = s.substr(0, p1);
        var base = s.substr(p1 + 1, p2 - p1 - 1);
        var linker = s.substr(p2 + 1);

        var ms = org.helm.webeditor.Monomers.getMonomer(org.helm.webeditor.HELM.SUGAR, org.helm.webeditor.IO.trimBracket(sugar));
        var ml = org.helm.webeditor.Monomers.getMonomer(org.helm.webeditor.HELM.LINKER, org.helm.webeditor.IO.trimBracket(linker));
        var mb = org.helm.webeditor.Monomers.getMonomer(org.helm.webeditor.HELM.BASE, org.helm.webeditor.IO.trimBracket(base));

        var m1 = org.helm.webeditor.Interface.createMol(org.helm.webeditor.monomers.getMolfile(ms));
        var m2 = org.helm.webeditor.Interface.createMol(org.helm.webeditor.monomers.getMolfile(ml));
        var m3 = org.helm.webeditor.Interface.createMol(org.helm.webeditor.monomers.getMolfile(mb));

        this.mergeMol(m1, "R2", m2, "R1");
        this.mergeMol(m1, "R3", m3, "R1");

        return m1.getMolfile();
    },

    /**
    * Cap R group (internal use)
    * @function capRGroup
    */
    capRGroup: function (m, r, mon) {
        var cap = mon == null || mon.at == null ? null : mon.at[r];
        if (cap == "OH")
            cap = "O";
        else if (cap != "H" && cap != "X" && cap != "O")
            return false;

        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.a1.alias == r || b.a2.alias == r) {
                var a = b.a1.alias == r ? b.a1 : b.a2;
                m.setAtomType(a, cap);
                a.alias = null;
                return true;
            }
        }
        return false;
    },

    /**
    * Find R group (internal use)
    * @function findR
    */
    findR: function (m, r1, a1) {
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.a1.alias == r1 && (a1 == null || b.a1._helmgroup == a1))
                return { b: b, a0: b.a2, a1: b.a1 };
            else if (b.a2.alias == r1 && (a1 == null || b.a2._helmgroup == a1))
                return { b: b, a0: b.a1, a1: b.a2 };
        }
        return null;
    },

    /**
    * Merge two molecule (internal use)
    * @function mergeMol
    */
    mergeMol: function (m, r1, src, r2, a1, a2) {
        this.joinMol(m, r1, src, r2, a1, a2);

        m.atoms = m.atoms.concat(src.atoms);
        m.bonds = m.bonds.concat(src.bonds);
        return m.getMolfile();
    },

    joinMol: function (m, r1, src, r2, a1, a2) {
        var t = this.findR(m, r1, a1);
        var s = this.findR(src, r2, a2);

        if (t != null && s != null) {
            this.extendDistance(t.a0.p, t.a1.p, 1);
            this.extendDistance(s.a0.p, s.a1.p, 1);

            // align
            src.offset(t.a1.p.x - s.a0.p.x, t.a1.p.y - s.a0.p.y);
            var deg = t.a1.p.angleAsOrigin(t.a0.p, s.a1.p);
            src.rotate(t.a1.p, -deg);

            // merge
            m.atoms.splice(scil.Utils.indexOf(m.atoms, t.a1), 1);
            src.atoms.splice(scil.Utils.indexOf(src.atoms, s.a1), 1);
            src.bonds.splice(scil.Utils.indexOf(src.bonds, s.b), 1);

            if (t.b.a1 == t.a1)
                t.b.a1 = s.a0;
            else
                t.b.a2 = s.a0;
        }
    },

    /**
    * Extend a point (internal use)
    * @function extendDistance
    */
    extendDistance: function (p0, p, s) {
        var dx = p.x - p0.x;
        var dy = p.y - p0.y;

        p.x = p0.x + s * dx;
        p.y = p0.y + s * dy;
    },

    /**
    * Create the popup dialog (internal use)
    * @function create
    */
    create: function () {
        if (this.dlg != null)
            return;

        var fields = {
            jsd: { type: "jsdraw", width: 180, height: 130, scale: this.molscale, viewonly: true },
            rs: { type: "html", viewonly: true, style: { borderTop: "solid 1px gray", width: 180} },
            notes: { type: "html", viewonly: true, style: { width: 180, color: "gray"} }
        };
        this.dlg = scil.Form.createDlgForm("", fields, null, { hidelabel: true, modal: false, noclose: true });
        this.dlg.hide(true);

        this.dlg.dialog.style.backgroundColor = "#fff";
        this.dlg.dialog.titleElement.style.borderBottom = "solid 1px #ddd";
        this.dlg.dialog.titleElement.style.textAlign = "center";
    },

    /**
    * Hide the dialog (internal use)
    * @function hide
    */
    hide: function () {
        this.clearTimer();
        if (this.dlg != null && this.dlg.isVisible()) {
            this.dlg.hide(true);
        }
    }
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* Formula class
* @class scil.helm.Formula
*/
scil.helm.Formula = {
    /**
    * Calculate the MF of a molecule (internal use)
    * @function getMF
    */
    getMF: function (m, html) {
        var stats = this.getAtomStats(m);
        if (stats == null)
            return null;

        var s = "";
        if (stats["C"] != null)
            s += "C" + this.subscription(stats["C"], html);
        if (stats["H"] != null)
            s += "H" + this.subscription(stats["H"], html);
        if (stats["N"] != null)
            s += "N" + this.subscription(stats["N"], html);
        if (stats["O"] != null)
            s += "O" + this.subscription(stats["O"], html);

        for (var e in stats) {
            if (e != "R" && e != "C" && e != "H" && e != "O" && e != "N")
                s += e + this.subscription(stats[e], html);
        }
        return s;
    },

    /**
    * Create subscription (internal use)
    * @function 
    */
    subscription: function (n, html) {
        if (n == 1)
            return "";
        return html ? "<sub>" + n + "</sub>" : n;
    },

    /**
    * Calculate the MW of a molecule (internal use)
    * @function getMW
    */
    getMW: function (m) {
        var stats = this.getAtomStats(m);
        if (stats == null)
            return null;

        var sum = 0;
        for (var e in stats) {
            if (e != "R")
                sum += stats[e] * scil.helm.Interface.getElementMass(e);
        }
        return Math.round(sum * 10000) / 10000.0;
    },

    /**
    * Calculate atom counts (internal use)
    * @function getAtomStats
    */
    getAtomStats: function (m) {
        var stats = {};
        this.getAtomStats2(m, stats);
        return stats;
    },

    getAtomStats2: function (m, stats) {
        var brackets = [];
        for (var i = 0; i < m.graphics.length; ++i) {
            var br = JSDraw2.Bracket.cast(m.graphics[i]);
            if (br != null && br.atoms != null && br.atoms.length > 0) {
                var n = br.getSubscript(m);
                if (!(n > 0))
                    return null;
                if (n > 1)
                    brackets.push({ br: br, n: n });
            }
        }

        var atoms = [];
        var list = [];
        for (var i = 0; i < m.atoms.length; ++i) {
            var a = m.atoms[i];
            if (a.elem == "?")
                return null;

            if (scil.helm.isHelmNode(a)) {
                if (a.biotype() == scil.helm.HELM.BLOB && a.superatom != null) {
                    //group
                    this.getAtomStats2(a.superatom, stats);
                }
                else {
                    list.push(a);
                    var n = this.getRepeat(brackets, a);
                    if (n > 0) {
                        for (var k = 1; k < n; ++k)
                            list.push(a);
                    }
                }
            }
            else {
                atoms.push(a);
            }
        }

        // chemistry
        var ret = atoms.length == null ? null : scil.helm.Interface.getAtomStats(m, atoms);
        JSDraw2.FormulaParser.mergeStats(stats, ret);

        if (list.length == 0)
            return stats;

        for (var i = 0; i < list.length; ++i)
            this.countMonomer(stats, scil.helm.Monomers.getMonomer(list[i]));

        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (scil.helm.isHelmNode(b.a1))
                this.deductR(stats, scil.helm.Monomers.getMonomer(b.a1), b.r1);
            if (scil.helm.isHelmNode(b.a2))
                this.deductR(stats, scil.helm.Monomers.getMonomer(b.a2), b.r2);
        }

        return stats;
    },

    getRepeat: function (brackets, a) {
        for (var i = 0; i < brackets.length; ++i) {
            var b = brackets[i];
            if (scil.Utils.indexOf(b.br.atoms, a) >= 0)
                return b.n;
        }

        return 1;
    },

    /**
    * Count monomers (internal use)
    * @function countMonomer
    */
    countMonomer: function (ret, m) {
        if (m == null)
            return;

        if (m.stats == null) {
            m.stats = scil.helm.Interface.molStats(scil.helm.monomers.getMolfile(m));
            for (var r in m.at) {
                var s = m.at[r];
                if (s == "H" || s == "OH") {
                    if (m.stats["H"] == null)
                        m.stats["H"] = 1;
                    else
                        ++m.stats["H"];
                }

                if (s == "OH") {
                    if (m.stats["O"] == null)
                        m.stats["O"] = 1;
                    else
                        ++m.stats["O"];
                }
            }
        }

        for (var e in m.stats) {
            if (ret[e] == null)
                ret[e] = m.stats[e];
            else
                ret[e] += m.stats[e];
        }
    },

    /**
    * Deduct R group (internal use)
    * @function deductR
    */
    deductR: function (ret, m, r) {
        if (m == null || m.at == null)
            return;

        var s = m.at["R" + r];
        if (s == "H") {
            --ret["H"];
        }
        else if (s == "OH") {
            --ret["H"];
            --ret["O"];
        }
    }
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* ExtinctionCoefficient class
* @class org.helm.webeditor.ExtinctionCoefficient
*/
org.helm.webeditor.ExtinctionCoefficient = {
    // εcalc = x(5500 M-1 cm-1) + y(1490 M-1 cm-1) + z(125 M-1 cm-1), where
    // “x” is the number of tryptophan residues per mole of protein, 
    // “y” is the number of tyrosine residues per mole of protein, 
    // "z” is the number of cystine residues per mole of protein.
    peptide: { W: 5500, Y: 1490, C: 62.5 },

    //Extinction Coefficient for nucleotide and dinucleotides in sngle strand E(260) M-1cm-1 x 10-3
    // Characterization of RNAs
    // [22] Absorbance Melting Curves of RNA, p304-325
    // Josehp Puglish and Ignacio Tinoco, Jr.
    // Methods in Enzymology, Vol. 180
    rna: {
        A: 15.34,
        C: 7.6,
        G: 12.16,
        U: 10.21,
        T: 8.7,
        AA: 13.65,
        AC: 10.67,
        AG: 12.79,
        AU: 12.14,
        AT: 11.42,
        CA: 10.67,
        CC: 7.52,
        CG: 9.39,
        CU: 8.37,
        CT: 7.66,
        GA: 12.92,
        GC: 9.19,
        GG: 11.43,
        GU: 10.96,
        GT: 10.22,
        UA: 12.52,
        UC: 8.90,
        UG: 10.40,
        UU: 10.11,
        UT: 9.45,
        TA: 11.78,
        TC: 8.15,
        TG: 9.70,
        TU: 9.45,
        TT: 8.61
    },

    /**
    * Calculate the extinction coefficient of a molecule (internal use)
    * @function calculate
    */
    calculate: function (m) {
        var chains = org.helm.webeditor.Chain.getChains(m);
        if (chains == null || chains.length == 0)
            return "";

        var sum = 0;
        for (var i = 0; i < chains.length; ++i) {
            var chain = chains[i];
            var ss = chain._getPolymers();
            for (var k = 0; k < ss.length; ++k) {
                if (ss[k].type == "RNA")
                    sum += this._calculateRNA(ss[k].atoms);
                else if (ss[k].type == "Peptide")
                    sum += this._calculatePeptide(ss[k].atoms);
            }
        }

        return sum;
    },

    /**
    * Calculate the extinction coefficient of a peptide (internal use)
    * @function _calculatePeptide
    */
    _calculatePeptide: function (atoms) {
        if (atoms == null || atoms.length == 0)
            return 0;

        var counts = {};
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            var m = org.helm.webeditor.Monomers.getMonomer(a);
            var e = m == null ? null : m.na;
            if (e != null && this.peptide[e]) {
                if (counts[e] == null)
                    counts[e] = 1;
                else
                    ++counts[e];
            }
        }

        var result = 0;
        for (var k in counts)
            result += this.peptide[k] * counts[k];
        return result / 1000.0;
    },

    /**
    * Calculate the extinction coefficient of a RNA (internal use)
    * @function _calculateRNA
    */
    _calculateRNA: function (atoms) {
        if (atoms == null || atoms.length == 0)
            return 0;

        var counts = {};
        var lastE = null;
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            var m = org.helm.webeditor.Monomers.getMonomer(a);
            var e = m == null ? null : m.na;
            if (e == null) {
                lastE = null;
                continue;
            }

            if (this.rna[e]) {
                if (counts[e] == null)
                    counts[e] = 1;
                else
                    ++counts[e];
            }
            else if (lastE != null && this.rna[lastE + e]) {
                if (counts[lastE + e] == null)
                    counts[lastE + e] = 1;
                else
                    ++counts[lastE + e];
            }

            lastE = e;
        }

        var result = 0;
        for (var k in counts)
            result += this.rna[k] * counts[k];
        return result;
    } 
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* HELM Editor App class
* @class org.helm.webeditor.App
*/
org.helm.webeditor.App = scil.extend(scil._base, {
    /**
    @property {MonomerExplorer} mex - Monomer Explorer
    **/
    /**
    @property {JSDraw2.Editor} canvas - Drawing Canvas
    **/
    /**
    @property {DIV} notation - HELM Notation
    **/
    /**
    @property {DIV} sequence - Biological Sequence
    **/
    /**
    @property {scil.Form} properties - HELM Property Table
    **/
    /**
    @property {JSDraw2.Editor} structureview - Structure Viewer
    **/

    /**
    * @constructor App
    * @param {DOM} parent - The parent element to host the Editor App
    * @param {dict} options - options on how to render the App
    * <pre>
    * mexfontsize: {string} Monomer Explorer font size, e.g. "90%"
    * mexrnapinontab: {bool} show RNA pin icon on its tab on Monomer Explorer
    * mexmonomerstab: {bool} show Monomers tab on Monomer Explorer
    * mexfavoritefirst: {bool} display favorite items first on Monomer Explorer
    * mexfilter: {bool} display Filter box on Monomer Explorer
    * sequenceviewonly: {bool} show Sequence View in viewonly mode
    * showabout: {bool} show about button
    * topmargin: {number} top margin
    * calculatorurl: {string} ajax web service url to calculate properties
    * cleanupurl: {string} ajax web service url to clean up structures
    * monomersurl: {string} ajax web service url to load all monomers
    * rulesurl: {string} ajax web service url to load all rules
    *
    * <b>Example:</b>
    *    &lt;div id="div1" style="margin: 5px; margin-top: 15px"&gt;&lt;/div&gt;
    *    &lt;script type="text/javascript"&gt;
    *     scil.ready(function () {
    *         var app = new scil.helm.App("div1", { showabout: false, mexfontsize: "90%", mexrnapinontab: true, 
    *             topmargin: 20, mexmonomerstab: true, sequenceviewonly: false, mexfavoritefirst: true, mexfilter: true });
    *     });
    *    &lt;/script&gt;
    * </pre>
    **/
    constructor: function (parent, options) {
        this.toolbarheight = 30;

        if (typeof (parent) == "string")
            parent = scil.byId(parent);
        this.mex = null;
        this.canvas = null;
        this.sequence = null;
        this.notation = null;
        this.properties = null;
        this.structureview = null;

        this.options = options == null ? {} : options;

        if (options.ambiguity != null)
            org.helm.webeditor.ambiguity = options.ambiguity;
        if (!scil.Utils.isNullOrEmpty(this.options.jsdrawservice))
            JSDrawServices = { url: this.options.jsdrawservice };

        if (this.options.monomercleanupurl != null && org.helm.webeditor.Monomers.cleanupurl == null)
            org.helm.webeditor.Monomers.cleanupurl = this.options.monomercleanupurl;

        if (this.options.rulesurl != null) {
            scil.Utils.ajax(this.options.rulesurl, function (ret) {
                if (ret.rules != null)
                    ret = ret.rules;
                org.helm.webeditor.RuleSet.loadDB(ret);
            });
        }

        if (this.options.monomersurl != null) {
            var me = this;
            scil.Utils.ajax(this.options.monomersurl, function (ret) {
                if (ret.monomers != null)
                    ret = ret.monomers;
                org.helm.webeditor.Monomers.loadDB(ret, me.options.monomerfun);
                me.init(parent);
            });
        }
        else {
            this.init(parent);
        }
    },

    /**
    * Calculate layout sizes (internal use)
    * @function calculateSizes
    */
    calculateSizes: function () {
        var d = dojo.window.getBox();
        if (this.options.topmargin > 0)
            d.h -= this.options.topmargin;

        var leftwidth = 0;
        if (this.page != null && this.page.explorer != null && this.page.explorer.left != null)
            leftwidth = this.page.explorer.left.offsetWidth;
        if (!(leftwidth > 0))
            leftwidth = 300;

        var ret = { height: 0, topheight: 0, bottomheight: 0, leftwidth: 0, rightwidth: 0 };
        ret.height = d.h - 90 - (this.options.mexfilter != false ? 30 : 0) - (this.options.mexfind ? 60 : 0);
        ret.leftwidth = leftwidth;
        ret.rightwidth = d.w - leftwidth - 25;
        ret.topheight = d.h * 0.7;
        ret.bottomheight = d.h - ret.topheight - 130;

        return ret;
    },

    /**
    * Intialize the App (internal use)
    * @function init
    */
    init: function (parent) {
        var me = this;
        var sizes = this.calculateSizes();

        var tree = {
            caption: this.options.topmargin > 0 ? null : "Palette",
            marginBottom: "2px",
            marginTop: this.options.topmargin > 0 ? "17px" : null,
            onresizetree: function (width) { me.resizeWindow(); },
            onrender: function (div) { me.treediv = div; me.createPalette(div, sizes.leftwidth - 10, sizes.height); }
        };
        this.page = new scil.Page(parent, tree, { resizable: true, leftwidth: sizes.leftwidth });
        scil.Utils.unselectable(this.page.explorer.left);

        var control = this.page.addDiv();
        var sel = scil.Utils.createSelect(control, ["Detailed Sequence", "Sequence"], "Detailed Sequence", null, { border: "none" });
        scil.connect(sel, "onchange", function () { me.swapCanvasSequence(); });

        this.canvasform = this.page.addForm({
            //caption: "Canvas",
            type: "custom",
            marginBottom: "2px",
            oncreate: function (div) { me.createCanvas(div, sizes.rightwidth, sizes.topheight); }
        });

        this.handle = this.page.addResizeHandle(function (delta) { return me.onresize(delta); }, 8);

        this.sequencebuttons = [
            { label: "Format", type: "select", items: ["", "RNA", "Peptide"], key: "format" },
            { src: scil.Utils.imgSrc("img/moveup.gif"), label: "Apply", title: "Apply Sequence", onclick: function () { me.updateCanvas("sequence", false); } },
            { src: scil.Utils.imgSrc("img/add.gif"), label: "Append", title: "Append Sequence", onclick: function () { me.updateCanvas("sequence", true); } }
        ];

        this.tabs = this.page.addTabs({ marginBottom: "2px", onShowTab: function () { me.updateProperties(); } });
        this.tabs.addForm({
            caption: "Sequence",
            type: "custom",
            tabkey: "sequence",
            buttons: this.options.sequenceviewonly ? null : this.sequencebuttons,
            oncreate: function (div) { me.createSequence(div, sizes.rightwidth, sizes.bottomheight); }
        });

        this.tabs.addForm({
            caption: "HELM",
            type: "custom",
            tabkey: "notation",
            buttons: this.options.sequenceviewonly ? null : [
                { src: scil.Utils.imgSrc("img/moveup.gif"), label: "Apply", title: "Apply HELM Notation", onclick: function () { me.updateCanvas("notation", false); } },
                { src: scil.Utils.imgSrc("img/add.gif"), label: "Append", title: "Append HELM Notation", onclick: function () { me.updateCanvas("notation", true); } },
                { src: scil.Utils.imgSrc("img/tick.gif"), label: "Validate", title: "Validate HELM Notation", onclick: function () { me.validateHelm(); } }
            ],
            oncreate: function (div) { me.createNotation(div, sizes.rightwidth, sizes.bottomheight); }
        });

        this.tabs.addForm({
            caption: "Properties",
            type: "custom",
            tabkey: "properties",
            oncreate: function (div) { me.createProperties(div, sizes.rightwidth, sizes.bottomheight); }
        });

        this.tabs.addForm({
            caption: "Structure View",
            type: "custom",
            tabkey: "structureview",
            oncreate: function (div) { me.createStructureView(div, sizes.rightwidth, sizes.bottomheight); }
        });

        scil.connect(window, "onresize", function (e) { me.resizeWindow(); });
    },

    validateHelm: function () {
        if (this.options.onValidateHelm != null) {
            this.options.onValidateHelm(this);
            return;
        }

        var url = this.options.validateurl;
        if (scil.Utils.isNullOrEmpty(url)) {
            scil.Utils.alert("The validation url is not configured yet");
            return;
        }

        this.setNotationBackgroundColor("white");
        var helm = scil.Utils.getInnerText(this.notation);
        if (scil.Utils.isNullOrEmpty(helm))
            return;

        var me = this;
        scil.Utils.ajax(url, function (ret) {
            me.setNotationBackgroundColor(ret.valid ? "#9fc" : "#fcf");
        }, { helm: helm });
    },

    /**
    * Resize Window (internal use)
    * @function resizeWindow
    */
    resizeWindow: function () {
        var sizes = this.calculateSizes();
        this.canvas.resize(sizes.rightwidth, sizes.topheight - 70);

        var s = { width: sizes.rightwidth + "px", height: sizes.bottomheight + "px" };
        scil.apply(this.sequence.style, s);
        scil.apply(this.notation.style, s);

        s = { width: sizes.rightwidth + "px", height: (sizes.bottomheight + this.toolbarheight) + "px" };
        scil.apply(this.properties.parent.style, s);

        this.structureview.resize(sizes.rightwidth, sizes.bottomheight + this.toolbarheight);

        this.mex.resize(sizes.height);
    },

    /**
    * Swap canvas and sequence view (internal use)
    * @function swapCanvasSequence
    */
    swapCanvasSequence: function () {
        var a = this.canvasform.form.dom;
        var h = this.handle;
        var b = this.tabs.tabs.dom;
        if (h.nextSibling == b) {
            a.parentNode.insertBefore(b, a);
            a.parentNode.insertBefore(h, a);
        }
        else {
            a.parentNode.insertBefore(b, a.nextSibling);
            a.parentNode.insertBefore(h, a.nextSibling);
        }
    },

    /**
    * Event handler when change window size (internal use)
    * @function onresize
    */
    onresize: function (delta) {
        if (this.handle.nextSibling == this.tabs.tabs.dom) {
            var top = this.canvas.dimension.y;
            var bottom = scil.Utils.parsePixel(this.sequence.style.height);
            if (top + delta > 80 && bottom - delta > 20) {
                this.canvas.resize(0, top + delta);
                this.sequence.style.height = (bottom - delta) + "px";
                this.notation.style.height = (bottom - delta) + "px";
                this.properties.parent.style.height = (bottom - delta) + "px";
                this.structureview.resize(0, bottom - delta);
                return true;
            }
        }
        else {
            var top = scil.Utils.parsePixel(this.sequence.style.height);
            var bottom = this.canvas.dimension.y;
            if (top + delta > 20 && bottom - delta > 80) {
                this.sequence.style.height = (top + delta) + "px";
                this.notation.style.height = (top + delta) + "px";
                this.properties.parent.style.height = (top + delta) + "px";
                this.structureview.resize(0, top + delta);
                this.canvas.resize(0, bottom - delta);
                return true;
            }
        }
        return false;
    },

    /**
    * create monomer explorer (internal use)
    * @function createPalette
    */
    createPalette: function (div, width, height) {
        var opt = scil.clone(this.options);
        opt.width = width;
        opt.height = height;
        this.mex = new org.helm.webeditor.MonomerExplorer(div, null, opt);
    },

    /**
    * create drawing canvas (internal use)
    * @function createCanvas
    */
    createCanvas: function (div, width, height) {
        div.style.border = "solid 1px #eee";

        var me = this;
        var args = {
            skin: "w8", showabout: this.options.showabout, showtoolbar: this.options.canvastoolbar != false, toolbarmode: "helm", showmonomerexplorer: true,
            inktools: false, width: width, height: height, ondatachange: function () { me.updateProperties(); },
            onselectionchanged: function () { me.onselectionchanged(); },
            onselectcurrent: function (e, obj, ed) { me.onselectcurrent(e, obj, ed); },
            onvalidatetext: function (s, editor) { return me.onvalidatetext(s, editor); }
        };

        this.canvas = org.helm.webeditor.Interface.createCanvas(div, args);
        this.canvas.helm.monomerexplorer = this.mex;
        this.mex.plugin = this.canvas.helm;

        this.canvas._testdeactivation = function (e, ed) {
            var src = e.target || e.srcElement;
            return scil.Utils.hasAnsestor(src, me.canvas.helm.monomerexplorer.div);
        };
    },

    onvalidatetext: function (s, editor) {
        if (scil.Utils.isNullOrEmpty(s))
            return;

        var t = editor.text;
        if (t != null && t.fieldtype == "BRACKET_TYPE" && t.anchors.length == 1 && JSDraw2.Bracket.cast(t.anchors[0]) != null) {
            if (!/^[*]|([0-9]+([-][0-9]+)?)$/.test(s)) {
                scil.Utils.alert("Invalid subscript");
                return false;
            }
        }
    },

    /**
    * Event handler when selecting an object (internal use)
    * @function onselectcurrent
    */
    onselectcurrent: function (e, obj, ed) {
        var a = JSDraw2.Atom.cast(obj);
        if (a == null || ed.start != null || ed.contextmenu != null && ed.contextmenu.isVisible()) {
            org.helm.webeditor.MolViewer.hide();
            return;
        }
        var type = a == null ? null : a.biotype();
        var set = org.helm.webeditor.Monomers.getMonomerSet(type);
        var s = a == null ? null : a.elem;
        var m = set == null ? null : set[scil.helm.symbolCase(s)];
        if (m != null && m.m == "" && a != null && a.superatom != null)
            m.m = a.superatom.getXml();

        org.helm.webeditor.MolViewer.show(e, type, m, s, ed, a);
    },

    /**
    * Create sequence view (internal use)
    * @function createSequence
    */
    createSequence: function (div, width, height) {
        var atts = {};
        if (!this.options.sequenceviewonly) {
            atts.contenteditable = "true";
            atts.spellcheck = "false";
        }
        this.sequence = scil.Utils.createElement(div, "div", null, { width: width, height: height, overfloatY: "scroll", wordBreak: "break-all" }, atts);
    },

    /**
    * create notation view (internal use)
    * @function createNotation
    */
    createNotation: function (div, width, height) {
        var atts = {};
        if (!this.options.sequenceviewonly) {
            atts.contenteditable = "true";
            atts.spellcheck = "false";
        }
        this.notation = scil.Utils.createElement(div, "div", null, { width: width, height: height, overfloatY: "scroll", wordBreak: "break-all" }, atts);
    },

    /**
    * Create property window (internal use)
    * @function createProperties
    */
    createProperties: function (div, width, height) {
        var d = scil.Utils.createElement(div, "div", null, { width: width, overflow: "scroll", height: height + this.toolbarheight });

        var fields = {
            mw: { label: "Molecular Weight" },
            mf: { label: "Molecular Formula" },
            ec: { label: "Extinction Coefficient" }
        };
        this.properties = new scil.Form({ viewonly: true });
        this.properties.render(d, fields, { immediately: true });
    },

    /**
    * Create structure view (internal use)
    * @function createStructureView
    */
    createStructureView: function (div, width, height) {
        var d = scil.Utils.createElement(div, "div", null, { width: width, height: height + this.toolbarheight });
        this.structureview = new JSDraw2.Editor(d, { viewonly: true })
    },

    /**
    * Resize Window (internal use)
    * @function resize
    */
    resize: function () {
        var d = dojo.window.getBox();
        var width = d.w;
        var height = d.h;
        var left = d.l;
        var top = d.t;
    },

    /**
    * Update Canvas from sequence/notation view (internal use)
    * @function updateCanvas
    */
    updateCanvas: function (key, append) {
        var format = null;

        var plugin = this.canvas.helm;
        var s = null;
        if (key == "sequence") {
            if (this.sequencebuttons != null)
                format = this.getValueByKey(this.sequencebuttons, "format");

            s = scil.Utils.trim(scil.Utils.getInnerText(this.sequence));
            if (/^((RNA)|(PEPTIDE)|(CHEM)|(BLOB))[0-9]+/.test(s)) {
                format = "HELM";
            }
            else {
                // fasta
                s = s.replace(/[\n][>|;].*[\r]?[\n]/ig, '').replace(/^[>|;].*[\r]?[\n]/i, '');
                // other space
                s = s.replace(/[ \t\r\n]+/g, '');
            }
        }
        else {
            s = scil.Utils.getInnerText(this.notation);
        }
        plugin.setSequence(s, format, plugin.getDefaultNodeType(org.helm.webeditor.HELM.SUGAR), plugin.getDefaultNodeType(org.helm.webeditor.HELM.LINKER), append);
    },

    /**
    * Tool function to get item by using its key (internal use)
    * @function getValueByKey
    */
    getValueByKey: function (list, key) {
        for (var i = 0; i < list.length; ++i) {
            if (list[i].key == key)
                return list[i].b.value;
        }
        return null;
    },

    /**
    * update helm properties (internal use)
    * @function updateProperties
    */
    updateProperties: function () {
        switch (this.tabs.tabs.currentTabKey()) {
            case "sequence":
                if (this.sequence != null)
                    this.sequence.innerHTML = this.canvas.getSequence(true);
                break;
            case "notation":
                if (this.notation != null) {
                    var helm = scil.Utils.getInnerText(this.notation);
                    var s = this.canvas.getHelm(true);
                    if (helm != s) {
                        this.notation.innerHTML = s;
                        this.setNotationBackgroundColor("white");
                    }
                }
                break;
            case "properties":
                this.calculateProperties();
                break;
            case "structureview":
                this.updateStructureView();
                break;
        }
    },

    setNotationBackgroundColor: function (c) {
        if (this.notation != null)
            this.notation.parentNode.parentNode.style.background = c;
    },

    /**
    * Event handler when selection is changed (internal use)
    * @function onselectionchanged
    */
    onselectionchanged: function () {
        switch (this.tabs.tabs.currentTabKey()) {
            case "sequence":
                if (this.sequence != null) {
                    this.sequence.innerHTML = this.canvas.getSequence(true);
                }
                break;
            case "notation":
                if (this.notation != null)
                    this.notation.innerHTML = this.canvas.getHelm(true);
                break;
            case "structureview":
                this.updateStructureView();
                break;
        }
    },

    /**
    * Calaulte helm structure properties (internal use)
    * @function calculateProperties
    */
    calculateProperties: function () {
        if (this.properties == null)
            return;

        var data = {};
        this.properties.setData(data);
        if (this.options.calculatorurl != null) {
            var me = this;
            var helm = this.canvas.getHelm();
            if (helm != null) {
                scil.Utils.ajax(this.options.calculatorurl, function (ret) {
                    me.properties.setData(ret);
                }, { helm: helm });
            }
        }
        else {
            data.mw = Math.round(this.canvas.getMolWeight() * 100) / 100;
            data.mf = this.canvas.getFormula(true);
            data.ec = Math.round(this.canvas.getExtinctionCoefficient(true) * 100) / 100;
            this.properties.setData(data);
        }
    },

    /**
    * Get selection as a molfile (internal use)
    * @function getSelectedAsMol
    */
    getSelectedAsMol: function (m) {
        var ret = new JSDraw2.Mol();
        for (var i = 0; i < m.atoms.length; ++i) {
            if (m.atoms[i].selected)
                ret.atoms.push(m.atoms[i]);
        }

        var atoms = ret.atoms;
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (b.selected && b.a1.selected && b.a2.selected)
                ret.bonds.push(b);
        }

        return ret;
    },

    /**
    * Tool function to select bonds of all selected atoms (internal use)
    * @function selectBondsOfSelectedAtoms
    */
    selectBondsOfSelectedAtoms: function (m) {
        var n = 0;
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            if (!b.selected && b.a1.selected && b.a2.selected) {
                b.selected = true;
                ++n;
            }
        }

        return n;
    },

    containsAll: function (list, sublist) {
        for (var i = 0; i < sublist.length; ++i) {
            if (scil.Utils.indexOf(list, sublist[i]) < 0)
                return false;
        }

        return true;
    },

    expandRepeat: function (br, repeat, m, selected) {
        var r1 = null;
        var r2 = null;
        var b1 = null;
        var b2 = null;

        var oldbonds = [];
        for (var i = 0; i < m.bonds.length; ++i) {
            var b = m.bonds[i];
            var i1 = scil.Utils.indexOf(br.atoms, b.a1);
            var i2 = scil.Utils.indexOf(br.atoms, b.a2);
            if (i1 >= 0 && i2 >= 0) {
                oldbonds.push(b);
            }
            else if (i1 >= 0 && i2 < 0) {
                if (b.r1 == 1) {
                    r1 = b.a1;
                    b1 = b;
                }
                else if (b.r1 == 2) {
                    r2 = b.a1;
                    b2 = b;
                }
            }
            else if (i2 >= 0 && i1 < 0) {
                if (b.r2 == 1) {
                    r1 = b.a2;
                    b1 = b;
                }
                else if (b.r2 == 2) {
                    r2 = b.a2;
                    b2 = b;
                }
            }
        }

        for (var count = 0; count < repeat - 1; ++count) {
            var newatoms = [];
            var na1 = null;
            var na2 = null;
            for (var i = 0; i < br.atoms.length; ++i) {
                var a0 = br.atoms[i];
                var na = a0.clone();
                selected.atoms.push(na);
                newatoms.push(na);

                if (a0 == r1)
                    na1 = na;
                else if (a0 == r2)
                    na2 = na;
            }

            for (var i = 0; i < oldbonds.length; ++i) {
                var b0 = oldbonds[i];
                var nb = b0.clone();
                selected.bonds.push(nb);
                m.bonds.push(nb);
                nb.a1 = newatoms[scil.Utils.indexOf(br.atoms, b0.a1)];
                nb.a2 = newatoms[scil.Utils.indexOf(br.atoms, b0.a2)];
            }

            var nb = null;
            if (b1 == null) {
                nb = new JSDraw2.Bond(null, null, JSDraw2.BONDTYPES.SINGLE);
                nb.r1 = 2;
                nb.r2 = 1;
            }
            else {
                nb = b1.clone();
            }
            selected.bonds.push(nb);
            m.bonds.push(nb);
            if (b1 == null || b1.a1 == r1) {
                nb.a1 = na1;
                nb.a2 = r2;
            }
            else {
                nb.a2 = na1;
                nb.a1 = r2;
            }

            if (b2 != null) {
                b2.replaceAtom(r2, na2);
                r2 = na2;
            }
        }
    },

    /**
    * Update structure view from Canvas (internal use)
    * @function updateStructureView
    */
    updateStructureView: function () {
        if (this.structureview == null)
            return;
        this.structureview.clear(true);

        var m2 = this.canvas.m.clone();
        for (var i = 0; i < m2.atoms.length; ++i)
            m2.atoms[i].__mol = null;

        if (this.selectBondsOfSelectedAtoms(m2) > 0)
            this.canvas.refresh();
        var selected = this.getSelectedAsMol(m2);
        if (selected == null || selected.atoms.length == 0)
            return;

        var atoms = selected.atoms;
        var bonds = selected.bonds;

        for (var i = 0; i < m2.graphics.length; ++i) {
            var br = JSDraw2.Bracket.cast(m2.graphics[i]);
            if (br == null)
                continue;

            var repeat = br == null ? null : br.getSubscript(m2);
            var n = parseInt(repeat);
            if (n > 1 && br.atoms != null && br.atoms.length > 0 && this.containsAll(atoms, br.atoms)) {
                this.expandRepeat(br, n, m2, selected)
            }
        }

        var bondlength = null;
        var mols = [];
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            var mon = org.helm.webeditor.Monomers.getMonomer(a);
            var m = org.helm.webeditor.Interface.createMol(org.helm.webeditor.monomers.getMolfile(mon));
            a.__mol = m;

            // cap R groups
            var connected = m2.getAllBonds(a);
            var used = {};
            for (var k = 0; k < connected.length; ++k) {
                var bond = connected[k];
                if (bond.a1 == a)
                    used["R" + bond.r1] = true;
                else if (bond.a2 == a)
                    used["R" + bond.r2] = true;
            }

            for (var r in mon.at) {
                if (used[r])
                    continue;

                this._replaceR(m, r, mon.at[r]);
            }

            if (m != null && !m.isEmpty()) {
                var d = m.medBondLength();
                if (!(bondlength > 0))
                    bondlength = d;
                else
                    m.scale(bondlength / d, new JSDraw2.Point(0, 0));
            }
            mols.push(m);
        }

        while (atoms.length > 0) {
            var a = atoms[0];
            atoms.splice(0, 1);
            this.connectNextMonomer(a, atoms, bonds);
        }

        var mol = mols[0];
        for (var i = 1; i < mols.length; ++i)
            mol.mergeMol(mols[i]);

        if (mol == null)
            return;

        if (this.options.cleanupurl != null) {
            if (this.options.onCleanUpStructure != null) {
                this.options.onCleanUpStructure(mol, this);
            }
            else {
                var me = this;
                scil.Utils.ajax(this.options.cleanupurl, function (ret) {
                    me.structureview.setMolfile(ret == null ? null : ret.output);
                }, { input: mol.getMolfile(), inputformat: "mol", outputformat: "mol" });
            }
        }
        else {
            this.structureview.setMolfile(mol.getMolfile());
        }
    },

    _replaceR: function (m, r, e) {
        if (e == "OH")
            e = "O";
        if (e != "H" && e != "O" && e != "X")
            return false;

        for (var i = 0; i < m.atoms.length; ++i) {
            var a = m.atoms[i];
            if (a.elem == "R" && a.alias == r) {
                a.elem = e;
                a.alias = null;
                return true;
            }
        }

        return false;
    },

    connectNextMonomer: function (a, atoms, bonds) {
        var m1 = a.__mol;
        var oas = [];
        for (var i = bonds.length - 1; i >= 0; --i) {
            var b = bonds[i];

            var r1 = null;
            var r2 = null;
            var oa = null;
            if (b.a1 == a) {
                r1 = b.r1 == null ? null : "R" + b.r1;
                r2 = b.r2 == null ? null : "R" + b.r2;
                oa = b.a2;
            }
            else if (b.a2 == a) {
                r1 = b.r2 == null ? null : "R" + b.r2;
                r2 = b.r1 == null ? null : "R" + b.r1;
                oa = b.a1;
            }
            if (oa == null)
                continue;

            bonds.splice(i, 1);
            if (oa.__mol == null)
                continue;

            var m2 = oa.__mol;
            if (r1 != null && r2 != null)
                org.helm.webeditor.MolViewer.joinMol(m1, r1, m2, r2);
            oas.push(oa);
        }

        for (var i = 0; i < oas.length; ++i) {
            var oa = oas[i];
            var p = scil.Utils.indexOf(atoms, oa);
            if (p == -1 || p == null)
                continue;

            atoms.splice(p, 1);
            this.connectNextMonomer(oa, atoms, bonds);
        }
    }
});

﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* AppToolbar class
* @class org.helm.webeditor.AppToolbar
*/
org.helm.webeditor.AppToolbar = scil.extend(scil._base, {
    constructor: function (parent, imgpath, buttons) {
        if (typeof(parent) == "string")
            parent = scil.byId(parent);

        this.div = scil.Utils.createElement(parent, "div", null, { position: "absolute", zIndex: 1, top: "-70px", width: "100%", height: 80, background: "#eee", borderBottom: "1px solid gray" });
        var tbody = scil.Utils.createTable(this.div, null, null, { width: "100%" });
        var tr = scil.Utils.createElement(tbody, "tr");

        scil.Utils.createElement(tr, "td", "<h2>" + scil.helm.kVersion + "</h2>", { width: "30%" });
        var td = scil.Utils.createElement(tr, "td", null, { width: "40%", textAlign: "center" });
        scil.Utils.createElement(tr, "td", null, { width: "30%" });

        tbody = scil.Utils.createTable(td, null, null, { textAlign: "center" });
        tbody.parentNode.setAttribute("align", "center");
        var tr1 = scil.Utils.createElement(tbody, "tr");
        var tr2 = scil.Utils.createElement(tbody, "tr");

        for (var i = 0; i < buttons.length; ++i) {
            var b = buttons[i];
            scil.Utils.createElement(tr2, "td", b.label, { padding: "0 10px 0 10px" });

            var d = scil.Utils.createElement(tr1, "td", null, { padding: "0 10px 0 10px" });
            if (b.url == null)
                d.innerHTML = org.helm.webeditor.AppToolbar.Resources.img(b.icon);
            else
                d.innerHTML = "<a href='" + b.url + "'>" + org.helm.webeditor.AppToolbar.Resources.img(b.icon) + "</a>";
        }

        var me = this;
        scil.connect(this.div, "onmouseout", function (e) {
            me.div.style.top = "-70px";
        });
        scil.connect(this.div, "onmouseover", function (e) {
            me.div.style.top = "0";
        });
    }
});


/**
* Resources class
* @class org.helm.webeditor.AppToolbar.Resources
*/
org.helm.webeditor.AppToolbar.Resources = {
    img: function (key) {
        var s = this[key];
        return "<img src='data:image/png;base64," + s + "'>";
    },

    'arrow.png': 'iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAIAAAAC64paAAAABGdBTUEAAK/INwWK6QAAABl0RVh0U29mdHdhcmUAQWRvYmUgSW1hZ2VSZWFkeXHJZTwAAAHvSURBVHjaYvz//z8DuQAggJjQ+H9+/fjw5tXLDz//EmEmQACha/78+unJ/Vu2H77x5MMvgvoBAogFjf/zx8/XD66ff/3q519GHycNST42JtyaAQIIi9R/hv8/P7+4emTTxr23X3/9g8d6gADCZe7/X1/fXDu0evXO2+9//sOlHyCA8Djq/58fH67vX7Zk6+2vfxmw6gcIICa8IfL/759P1/bMm7vh5odfWOwHCCAmwhHy/8fN/YuXbrn1+tuff6gyAAHERExi+Pf3y43DK9bsvPniM4oDAAKIibi0BPL/jaPrN+29iRz/AAHERHRaBIX/9WNbdp+49/Lzb4h2gABiIUYfM6eMtr4sFzMjIwOLiBAXOzMjRBwggFiIczSbgr6Nrjg7KxMjKxcvHwdUN0AA4XI2UJGgkpGtgSIrKMD+vH76+j+/kLC4uIgwLzsrzGaAAMKimZGBkZ1HRNXIydXJ0tJQlQOo8v/nO+evP/vyGy2rAAQQumY2dnYhaVUtY1t7W2NdBTFFTX1lQZCaT0/P33j4+QdqTAMEELqfeQTFdM0dlDmkVKS5mRn+cYmqGOuI3Tj0ioOX7T/QXlSbAQKIEX9J8v/fnw/3T2w7/opXUklLT0NeGBRmcFmAAGIkXAz9//X58z8uXlgQIwGAAAMAMJ7FYOAyZgwAAAAASUVORK5CYII=',
    'canvas-1.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAAFx0lEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAATTA1gMEEE7rWX+/Zf9yjXIL2H4+Z/9yHZcsQADhtv7DRa6rHcz/f1JiNyPDf9aX+9kercWlACCAWHBq/feL4c9Xhv9/gWYw//vO+hnog/8kh+3XJ2x35vwVtcClACCAcFr/n5mD4d8fxn+/GZgYmH++5LzSxvCfZOsZfn9i+Pr0vwQXLnmAAMJtPSsf0OtMf78xsPD/5pT9YjoF4XlGiApULlwEicv66Sr7xcb/HOK4bAEIIJzW/2MTAprE9PM1A7vkfwbm32wSZMT9f0F2djb+f7xKuBQABBBe65k5mL/cYeDTAydRcEiQmuzfnGD4++MPjwouBQABhNP6v8y8vyXdmN5dYpQK+M/AxP7uFMf1XmSPIYUyDu7fnwxfH/8TtfjDIYXLFoAAwh33DIx/RMw5z5UKPF3GIO/H8OcVg7ABqgV4ox2YZd5dZfjz/Y+Y9T9GVly2AAQQC75ky6PKJqjP/HQvg7Qzg1oc2FhYemNkxMcFgq9PGa5M/M/I8kvCFY8VAAGEz/q/TJw/leK5TmYy3JzHoFfCwCEKtuMvw4/XoFLhP1LSZ2RA4QLVPNjI8PzAb/Us/GkWIIBY8KedX7xaLJp5bNf6GFi4GbSyGDiEQWXR5QkMn+/jjDSg3b8+M3x59Ecu8LukD37zAQKIgPXAFPBNOojx53vWG1MYfrxhMKhgYGJh+PIEFLZYcgvQ4vcM354z/P/zVyXhm3r+P0Y2/OYDBBAjka0drueb2U8VMAjrMpi0MYgagxzByIQU32BPf7jOcLGb4cmeP+qpX1Xz8KQ4OAAIIEbiG1tsPx5z7/Vk+P2FQT2JwbCSgUMEEd/AwvX6LIbrMxh+f/ths/C7oDmRZgIEECNJbT1gouJ+uJD1bCUDCw+DTjaDRgbIEUzMDDfmAP39R87vm0bpX0Z24g0ECCBGMpqazP9/AOtilltzGNj4GBT8GNQSGHiVPvxk/Y9S8hAFAAKIkeyWLvO/b3yXqhiApcL/fwySNp8M+/8ycZFqCEAAkdzY4nxzQOCgJ/v3+yDLzNoZ3NYyqMUyvLvCt9uB++kaYNVAkmkAAUSC71n+fOS60sj8/ACDhDWDYQ2DgBoo6UFKoU/3GO6tZni27x+v8ne9xl+c8kSaCRBAxFrP8ekKx4Vaxn/fGTRSGCTtGHgVGRiZESXuvz+gkuDVCYbbS/7/ePtLu+S7qBMxSQEggIiynvPtIY7ThaDyBFhzixgysPFjaxzyMSgGM/z7yXBnKcOrU79UU75LB/9jJFCsAQQQC2G735/kOF/FwMrLICLJwMIJCmpg8YfS7gOV+KAagZmFQUSfgV2QgW0J290FwOABlpjAyhqP4QABxEKoqHnKdnMKg5Aeg3IkAwsHSp2GVssBXcYtC4oRfnUG9QRgmmC7M/cflyz+IggggPBZz/T/D9uTdUx/fzBoljBI2oJLHdw1LLIIrzKDSiTDtxdsdxf9NtT6w8yLywqAAMJnPev3h6x3FzHwazB8f8HwYBP2Kglr+uKWYuBRYJB2Ybo+g+31kT8SnrisAAggfNazvDsDimaWRwzX55BWOEjaMCiFM4gYMPCrsbzczyjhgSsXAAQQC55Cjfn9BQYBTQZJe5JLRGAGASYFdiEGIR3muytYfr3E1egACCDc1v96z/z+KoO8D4N+Bc7kRpALTowsn+/+FsZuPUAA4bSe6ddbhr/fGYAZ98tDcqoEVh5Q8QDMhKx8zF9uMwhbY1UFEEA4rWf8/ZHh2zOG5wcZvj7B0a/BywWWTtKuIEewcjECyyscACCAcFv/5xuoEQFswODWjK+RyK/77TcHA5uqADMX45/PuJQBBBDeYkfM4rP5PDy5lojmyb/fgkZMXx/gUgAQQLh7Odxyv+TD/zJzU9K/B5a4fwV0/rMJ4lIAEECMAzuwBhBAAzy2AxBAA2w9QIABAKhg/GFIL73sAAAAAElFTkSuQmCC',
    'canvas-2.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAADt0lEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAATTA1gME0ABbDxBAA2w9QAANsPUAATTA1gMEEFHWP3j0jFRzf/789e79R4LKAAKIhRiznj59riAnBWT8+fOXGPVAZY+fglwsJMiPXyVAABFl/X8gBIPb9x7ChMAijIxY2IyMf/78efvmrZSUJEGTAQKIKOt5uLn+/fvPxMSoqaZEZOC/fPUWqIWgMoAAIiruBQQEfv76Bfc2QfTz1+8fP3/x8/MQNBkggIjyvYSYyOu372WlxX/8+AWJVHwx9Z/hxctXoiLCXJyEAx8ggIiynoOD7cuXLwwM4kCGqrICfsV///779v0nPz8fMSYDBBCx+V5cTOTh4+fEJfs/IsIC4qLCxCgGCCAWIq0HZqH3Hz4BXSAvCwrSL1+/vf/w+d+/f2gW//r9G8iQl5FiZibKYwABREKpp6woC7Tg8dOXQDYnBwcLCwsbDDAzM3/79uPps+cfP36SlZbg4uIg0kyAAGIktbH17MVrRgZGQQE+BkYmIMnGxgosZJ69ePXq1Rtubk5hIQF2NlbiTQMIIEYy2nrA9P/x0+eXr9/x8/EKCwsCS4U/v38Bcxo3Fyek7CEeAAQQI3lNzf///9++++jjp0/MzCyCAvzyshJMTOTUXgABxEKW3aBSXVVFEZgAnz9/+ebtOx5uTmBEsLAwk2oUQACR7HtgEnvx6s3Pn7811JX5+XiAJevbdx+ePX8JtJubiwPoCEZSIgAggEiz/vmLNxcuXZaXk5eSFOfk5GBkgtrExsIEjI43bz98/vJVXpbYXAcEAAFEgvVPnr28ceMWDy8vBwc7mgUc7OyQ0vDDx8+vXr9VlJdhZSUqWgECiFjrgdF8/+ETbQ1VJiZGtHSAFtjA2Hn15h2wdCImFgACiNhQunPvobaGCprd8FoeGQDLHGChCywiiTEWIICICiJgBfrkyVMWZmAxR6BIkRAX5ePlBpaJX79+J8ZkgAAiyvoXL9+oq6nISBOuQCFFHjBlAIthYHOPnZ0Nv3qAACLK+k+fP2lrqGGGPD5zWZg/f/lG0HqAACLK+sePn3JzcRNsN4KSEhMTsAhiAoMvXz8DEwF+9QABRJT1f//+BRqHVr1iBf9AANgqBLqDEaiLoHqAACLKeilJCUhDm3jAysICLJcIKgMIIMaBHd0ACKAB7mQBBNAAWw8QQANsPUAADbD1AAHEQjujmf/9AHf/mP4yAwsf7EUWQADR0Pq/TIQzHkAADXDgAwQYAMXQUECyhsOVAAAAAElFTkSuQmCC',
    'helm20.png': 'iVBORw0KGgoAAAANSUhEUgAAAHYAAAAlCAIAAAAY34ofAAAACXBIWXMAAA7EAAAOxAGVKw4bAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAAKIElEQVR42mJ89+4dwyigJQAIIKbRIKA1AAig0SCmOQAIoNEgpjkACKDRIKY5AAig0SCmOQAIoNEgpjkACKDRIKY5AAig0SCmOQAIoNEgpjkACKDRIKY5AAig0SCmOQAIoNEgpjkACEAplSQBDIMgFf//4gj10MkpS+sNZRAd8sqQaBbu62EVAUzMGt3RgumZGfGqiCWLCbeLqcAxBOIoZeLLqS1XpHsA/4LVRmS+e8G1HgFEIIj/fHx2aM+Bd3wang5G3Kzosk9Oblm69655SIqDGg+Q++vFtdXrtt5+9v7n3/9glzEwAj0EDEhGIItD2dghwMNWkIPh77c3O9ZseMujFhLgyMX0F3vm+v/n6rFdO47fNXT0sjJUZsMaJv9+nt+5aseFTy7xaaZSrPjC6M/Pp/dunD59+trdp/+5RC2tbU0NNPk4CSevf7+/379x6dLVu9+Z+LR1NWQkJQT4uJmZSAtrgAAiYM37x/dvXrr0hI/RyEBLSZQDJVn+/3Xt4PH3X1j2H7vnoKYHFLm5f9e9Z29//GNlZ2dBdsf///+BYf3y1etP3/8IcrB8e3Uf6Nf//N/uvrPXFcFh8c93V65eev3mw84tWxiYA+z0FZgY/6OlydtHNqw/dPXXf4atp56aBijgDN9fX0/v3bjz8KWvwBQMBN8f7d689vxl4+QEP0G2//jC98eHozs37Ttx7TtY1aVzRxhZeZR1bUL87AU4SAhlgAAiEMT//oEQw38c+QecWhlAJQkDI8P/r59+AxUz8Um7OJiIC/EB0zAwAf8D6QeVEwIiYmKCIOv+M4LM+w8sZP7h9CEjEwsw/YMUf3++Z/tOPp4gIxWU2Hh549Ca7ed/QQ3A6WGgVY9vnD109PJ3Zk5eIUEBHlbGnx/fvPrw9umVDUeNEh2lcdZR//9cO3346OlrP1k4uDk5OVgY/vz6+ePH1zuXj61n44sPMGZi+E9kEAMEEAvB8gvsUkYGLCXRf5gHGcHEvy+fvgGFmGVNrc2N8WQmoAww8ICx8w9kAA51TCzMjMywwur+9u27uEL8NSS5ISJfnl1dvWLnR7gf8ZSSvz/fuHj83R9mPgVtKw9fZ02ur/dOblq47syn32/vPGPAHcTfX905e/b0h7+sQI0mpoZawv/ePr5z49L5S4++PLx58tY7Yw0hYlMxQAARKPjBJSkwLIFR+AtN6u/f37/+Iiei/7//Q2IWT+qEKGcGhwkwD+BOxYzAOgmUjEFlORPjl2fXdu7c+/TTH1CJ//HZ+oVLnvwCxgLhau7vr+8/f/5gZheWlFXVUhXjYgSlxV+g2ouRlYMVT9p/dPv+y9e/GDnlNPUt7G1NdA3M7ZydHK20hZj+/fz55fzNl8QXFAABRKjIZwSF37/3t9Yvmy3Ex/XzBzCzfPvx4+fv33/+gurn/wxIcfQPkujuH9246SUfN/t/EIBEFKuorIqOmjwb039wIcAEDJrfoEjCHRWgehKUd3jEpfiZ/n579eL5zXPbtzE7WOte2zTv6kcGFk5uMXUVzqsX7v5GZCcsSYRbzMozWv49o7yKIt/vj68enjuwbc/lTwwsPFxqJqo4df35evPOjQ9/GTllZaUVZETZQCmCiYOHT1RKnpPx9c/fH5+8ZmAQIzKIAQKIQBAzMoLD+O+vj29ffXyLO+cjJ+YvL0+fwIhkjjt/I0NM1YWQNeBJ7KBGEjh+2cS1zbQFft86svvM0zsXDgMRSD+7gLC6U7gX9+FrFwgYxMAoJK0iJA2tvvZs3HL00Q9Gdn51Mx9PTW7cNeSX7z9//GNgFhbkF+bngaYiRjZ2Th4BQUaGp39+fnj7F9j4+/+PmCAGCCAiUjE4CplZ2Xh5eYCAm5uLA9hiYGH89undo7sPvv7HDGrMVMGhrKOuIicIq8pAqRgYin9ABT0zntgFG8iko2nIqa7C8H/JtjMP/4JTppKmVVy4OQfDN2ik/iNc8/z99nbn0llHHn5nYBPQsfCOcNNgwF1f/QO2tkE5lIGVhYUN6FV4HczEwgpqP/7/8uULPqejAoAAYiGiLGZkFFCPjg7UkuVDbWx+29LZfPQruhZWncC6KDNW3PUdqHgFtyn+/cOXfVA4bLzmHsHv3sw69uCnuKKWl48dOyOo4GeCVcoEWvdf3+xbu+TY/U/A8NW1dA/z0MPfHmBi5WRlZWVk+Pb124+v33/+Z+AEu+bP79/fv34GdcR4ePiI770ABBDBVMwILo6BVQsLZuizYA0VJlZmUC7Hl3mZGMGNtv//Sag0eMS8U0ptPn5l5ebn4UDNN3jN+fvt3ZEtq4/efPmXmV/P0jXI3YiFkYC9wBhVlFO6dvfMu+cvnz9/801WlhsYJ7++fX738tHHf8C0zSUuRmQpAQQAAUSoLAaH3H9QF5JwcDDhKSvQajJwQYGvuvsPrSwZkZIzEwu7oDA7lkoAT8vk7/dLx3YevvL4F4uQrLZVkJsxOxMRpQojm7ycjCjPuYevb1w+JygkwG4iz/X98Y0Lpy89/c3IysupriNHfOIACCAWolT9J5BSUNz36OShY7/EhPiZIC0+SOgws/IJiYoK8jBC0j9Q2c/PN0/sUvinDml9Q5sebPxKCpIsTFA7GYiLMzzgx/sXt65c+/b7HyP7729vbl84x8HGwghqEQIBCws7l6CivASk0Pj06tlXJj5JER6IRhFlLW3t6y9O3np18/iOFzeOcTH///Hp/bsPv5m45OUNzBS5YA0owgAggAilYmZwKDFi9SwjIxO8XQdqA7BxsAHb7P8+Pd236yWoA42cAhmZ2ETUw0L8lETZmZhZgOXc/z8/n988u+LeeUiShVnH7RCRbKciCDOWEQIY8HelcJf7zMAylRWcbn9//fDi/vat9yFWgU1lYmLnMw/NcFPmYPr6eMfmLS9/CYXGhkjwgCtjVj4TW6e3bz+evf3i/Y+v74GuAY0DsInIqft4W3Mw/iM+mgECiLm8vByPNBsby9ev33hlVA20FDhYGNHGzpj/vrv14o+du7OCMDBwGfl5/9+5/+Dbz79A8AcGfkPAX0ZhGWUDHTUedkYWNs4/7x7cf/nhz1+oLEzt3/9swuaWpqI8rP8ZWblY/j1//cPa3lZWhAt7EDKyMP14ducds7OXsxw/C462DK+ktDTjr2/v3r3/+Qto0V8IANsGtE7I2sZClAvUBn947fInNkljfRV2ZqhtrFwCShoakkI8H9+8+PHrPwuvuK2Lh6+HnTgfG0k5CSAAp+a2AyAIguFhoL3/28ZJqWZeZHPrBTgNGP83YPnTdk0xzMFjdVFHGrUXPj60HSARdopYTSPh+r7S0lYyDeu2hJnG4jnTgrSEElUx896/cM8NJNwLPRjWT9k+cafCccRFbOHox+JqAohxdNkgrQFAAI3OetAcAATQaBDTHAAE0GgQ0xwABNBoENMcAATQaBDTHAAE0GgQ0xwABNBoENMcAAQYAIKsCm3IdpQWAAAAAElFTkSuQmCC',
    'monomers-1.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAAIBElEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAATTA1gME0ABbDxBAA2w9QAANsPUAATTA1gMEEGHrmf5/Zv73llRzmf+9Z/77ipHhP35lAAHEhNfiryyfN7M+CGV+EA426zcTwy8iLP4A0vu6j/lhFMuHZSw/rzEy/MGlGCCAWLAb8f8j44d1zF/2MXEa/OfUYfj9nP3zpv9/3v35/5dRKPkvsxhO497PY/xx9a9k7382RaZft1lZ2f+87mBklfkvEPeHXQNTPUAAsWALkN+MH9Yyf1zHwMTCIuLGxCH3/9/vvx+O/v16iunPcwYWGQZ+j79Momi6GBn+MX87wvjzJgMzP8uvu38Ek7jEvIF6/345z/Tj0r+3M5jFy/4yS6HpAgggxs/XMn5LtEE4rG/7GX7cA7vgJ6toNIugHSOrAAPDv7/f7jIycTCyCv/9eOb3qxV//34BqeZ1/sPn+5+BmfnvS6YXLSBdv++ziGezivoysvIw/PvHwMTE8A8Y7P//fDjx+9Xyv///MPz9zsDE9k+8+i+TENAKpi/7AAKI8fPlaAZmrv9CyYyftjL8esj44wojMye7ymxmLjWwk/79//vj1+PpDAx/WMRDmTmV/nw4/OvpxP+/HgHl/nMaMbAIMPx6zvjzOiMTOzN/ALtCAdACoIMZmJgZ/v0Ghh/EHf9/f/55v+nv5wNATf+5LRkYWRn+fWL6/RAggBg/vb7K9DAWGhSMzEz8XqziEUycCoxAFYxMQCGG/3///Xz569m8f19PswgFsoh4MPz7+fvlxr8f1vz/+xWki4mdiduKVSqZmUcT7OC/YMOALmCFhcFfoKJ/P5/+frHu78et//+A8hEjM/dfmWkAAcQIrO9Zf99gfJwDFOLUXs/EJgEKe0ZmoK0g60E5B+iCf////wKG5L9vN3+/2czw9wuLiP//f79+PawARgoTnxuHYgU41bDD/A3MIMxgESZwqoA5COj1/39/PZr2591yBon6X1y2AAHEhJwAGZn5fr1Y/vPhJKBLYXYDyX/AMAB6kZGFm5lXj0XQ6f/fb38/nGCA5mlGkFuByqBhzgQOc1ao3cC4h9gNNA0oyMjCyMgOJBlgmgECCGIxLF/+//v/x8v/f96AZf+DlfwBu4MRIg20iZnfgolTBeimf9/vIVI9IziOQSH/Bx7fIDYoFP+DTfjL8I8REsFgLiiigSRAAKEWOyChfyD0/x/Y35CQZ4A6H8QGhgMLE7sEE7s0WDEMANUzgsMJqAtiN8jFTCC7QYEHNgdk8X+Yt4EyoFABCCAWlLIPmDdAdvwHKQW6lxFsN4j9D5gDIP6Eeug/cmn6H5JCoVKQtAYJMFDiZQI5DqT8H9gosFthvgcIIFTrgUKMjFCzGGBhBfITI5jHCLIJEsjQ8ICDf9C0wghN59jCHOyx/3AngwBAAKEGPlT6PwM0kBhhYQ6WgsQFIyQM/sMSBAMsgsDplBGWSkAOAtJgWyF2QQ2EGw0CAAHEAvUIlxEDEzc4xULMZYT6G+IgRkao14Ex9v3e/5/PmLg1IWGIcAEDOL7+M0JjHeTXf7CAYYQZ+w/mJSgACCCQ9b9ZlTmkUhn+/QAmYBYRV6ANjGzCUKVQu8G++/f7349Hf97s/v/3PSurKAOiMv2PiKP/DNAwR+hlhAYn4z+on5GSDUAAsYAj4Nffdwf////OxC7DxKXGCMzBkAKLEar238/nfz+e/P/nIyO7NKhoAxYvjMz/vj8Am/X7/88nQAVMbOJg9f9ACegfLNFAEgHE38Bk8/v9vy8X/32/BrceIIBA1jP+//bn3QqQUWJh/z6eBJrOzKsPjA9ITP/9evPP231AHzNzazJxKjGyCP7/8fDP2z3/vp4Ap5Y/DD+v/n62kInPjIXPFFTZ/P8HtQ+S1mDB8/fzxT/vj/37eur/zzvAsvIfsxBQAiCAQNb/ZRRgEq9i/Hn3029xts/7mL6fZ/gfy8xnBKx7wGUlsJpiZuGxYuY3B8bO30/ASm/N/+9ngX5hZFdkZBb9//PK3w8b//24/u/7fXbpOFAZCk2AkDwMSjH/fjz+/WbLv0/AiPvGyGP3jz/kD5sqUAVAADHC+3jAav4fAyvr99NM7+YysQFjQYuBmYOZz4yJVfD/v+//fjwFGgEMin9fgaF3+b9AxH9O/X8sYv8ZOZiAFe6vxwwftzL+uMgingqszZjYZVkELMHFwL/fb3cAa1xgEfn/+/V/PA7/2TX+sSn+ZeSFWAoQQIyYXUyWvy+YXjYz/LgFYgsGMQvaAwv832+3/fu44/+/nyC3C8X8YVH4z8CE2kB6z/zlOMPrPlBYy878yyrPw/QIWHP+uJ3378cNRlA13/iH0xhNF0AAMWLt4bL8uc/8YTXDl4NA5zMCG1uMLP+BzZi/P/5zGv4TzfnDLImrscX2fjbjp80/5TexfVrL9O0Mi7Dn71dL/jOL/Ofz/M1li6keIIAYcXWwmRm+Mf24yfj1EOOXPcBoBur/z2X5j0MNmFDwtDOBKZXl9+3frBrsj8IZ/n78J17BwMT/j10e2ETEqh4ggBjx9++Z/39mfV7x/8/TP7IL/zLyE9/QZv1+iuHfl7/ctsD0hEcZQACx4DflP7B6ZuJkYORiYCCtQ/Kb04wYZQABRNDQ/zTt5QAEEJF+YqSR9QABxES0xTQJBoAAIhT3DEz/uKwYft8HNY1pAAACiHFgR7YAAmiAO9gAATTA1gME0ABbDxBgAHggOaqwTElCAAAAAElFTkSuQmCC',
    'monomers-2.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAADiElEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAATTA1gME0ABbDxBAA2w9QAANsPUAATTA1gME0ABbDxBAA2w9QAANsPUAATTA1gMEEGnW////nxg1//79I9JAgAAiwfq/f/8+efIEWEPid8SDBw/u3LlDpJkAAYTF+t+/f6OJfPjw4fLly3/+/GFjY/v58+f79++vX7/+69cvNF1XrlwBOk5RUVFYWBio5cKFC1++fMFvPUAAsqdQBYAQiqE8UTC8aFDB5v//h9XgF4hiMQjGGxeu3No2tvdGf6m1xswhBCHEvRcUu621RBRjPOfMOffepRTnXEpJSonFvXellNaaX6y1UFVrRTDnDAuvjDG898aY79YjgBixNje+f/9+6dKlb9++AcNZQEBAS0tLRESEmZkZ6Kf79+8DFQgJCQHDAMgGigCtB3I1NTWB/gaqgRsCTAGPHj0ChhPQxUBzDAwMxMTE0CwCCCBGPK2d3bt3q6qqKikpAX0AtAMeyM+ePQMmAj4wuHbtGtAaZ2dnZD8hA2BkHTp0CBgdysrKmLIAAYQv6QG9ArQYGAxAj8IFgcEoLy9vZGQE9ArQKZycnMAEgctuIABGB1AB3PVoACCAmIhJ8JhJHWgrMOKBYU5hvgcIIMLWAxMg5YUBLkMAAoiJ8qKGoPvwAIAAYsEjZ2dnB0nJaFELTE1Pnz4F5klgBiO+KMQqDhBA+KwH5hkuLi5gigWSEBFgyQNM9sAsDhSRkJAA5ihgxgP6HpjLgTkTqyHAnAVUxsvLi1UWIAAfdXACAAzCUHQO919KhEzgEn2nnkrvSowk/yfPHA0COINC4ra7TIOPsiWBB++R7e52ELxoyl03OTOIpCxV9ZQ4AoiRYC8HGM7A0grIAAYDsAgCGgr0PTDf8/DwAN0E8TTQ9/fu3QOGjYqKCrCoAIoAS6Tbt28DGUA1UlJSuNIHQAAxEtPJAvoSaBzQx0ALgHkdmI+BtgI9hJybgbJAZwEDDJhQgDENTB9ANwGLWKDX8ZgMEECMxPfxIF4EBgCw2AEGOJ4aD2g9sAAmxkyAAGIkqYsJNBrob2CwU6u5ARBApDU3gBbj8TcZACCASLOe+GYMkQAggFhIUg0s5IG1PhWtBwgg0nxPXbuBACCABrilCxBAA2w9QAANsPUAATTA1gME0ABbDxBAA2w9QAANsPUAATTA1gMEGAAj0m/VIoHc6gAAAABJRU5ErkJggg==',
    'settings-1.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAAHmklEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAATTA1gMEEAHrGRn+0dR6gADCZz3Ln+eMB/0Zv36hnfUAAcSE29//GZ8u/PfxOcPlObSzHiCAcFrP/OfJn9urGBj+MrxdzviFVgEAEEBMuKKc8flahl/szIIxjNz/GC7OpJH1AAGE3Xrm30/+3lrJyO3FrpvKph76/z0wAL7SwnqAAGLC7vVnq///EmQWCmSTlWCV9mbiZWC4MJ0W1gMEELr1jP//ML87/vfuZiZuazYlQ0YWNmZOdTaVoP/vVzC+fc74n8r5ECCAWJi/PmX4/Oj/m0MM7w4w/P74/9/fv3/+MjBKsYh5sYpzA93DyMrPKuv369baf8cDGNkYGZlYGJl4GET8GUTNGXik/nOJ/GdmIdt6gABifL9KleG3CCOHGhOXLCOnCBO3MBO3EBOXFLOgCjMfD1jN////vv15c/nf57f/vr399/XNv6+v//948u/zeYb/vxkUm/6ru5FtPUAAsTDxCP9985OZ25RN1ZVFRAgY2ozMbAwsrIxMzHAnMjJxs4qa/Bf8/f/v7/9/fjL8+fjz3tJfX44xsFoxylv/pyDwAQKIidFkAbOEzL8ve/68vPn/DwcTFx8jOwcjMzMwBaIGEwsjKycTBx8zN++fdyf/vNjBwBPAaN3+j4ObkrgHCCBGYHOD5c8rhiuN/959Z5GMZ1d3YhHmw6n8//ef99f8vDHlP4MNg1HRP3Yu3IXmX6ZPD0A1BiPzfw4xXCoBAogR0tph/vOa8Xr3vzdvWMSj2dVcWUT5sdr9496qX9en/GO0YTAp+8/GgS9Sf9z+eyiR4S8wDFkZxUL+GWdhVQYQQNCM95dFlEE9n5Hrw5/nu36/fMaAJT6BCfD1r9sz/r7nYjAuxbQbWEcw/f7E+O8PlPv95v8fX/7/+vz/1/v/LzfCDWT+epfxHyL3AgQQIt//ZZdg4pH4/+8DAzDj4aqFfr1kFHD+z86JGdTMb3YzHI9lvDid6csr5h/P/99fBrPy//+/HxjfP2X6+4P10ax/x6IYL0xj/AN1JUAAsSAVOMBU/ZWRXYqRkx+W7v4zAF3KCEmGjIyMwkz8gn8eHWBgyEUvpN/u+3ep5t/nXwzvFzK+2MrAwfHv+xOE9L9fDMdTGPj4fn+8CazC/n+bx8TA9t8wFWgkQAAhfM/099P/nx8YWcWYuAWAOv7//vTnzaUfd7b8vH/17+dvIK8wcjLzaTEwPGP8+w81/v4wvDn0//tvMO/P/18v/n168P/3H5SI+/P83zuQ3RDn/Hu2mfEvSAFAACEVWL9e/f/xk4lbhJH1z5931/68Pvv76f6/r64xsumySDuxSpmxCMsw8ekyMh1n+P2LgRkR9/8YWJjkYpheH/j77jMDUYUAG6Nh+38WViALIICQAv/Hs/8/GRnYv/15vfPv2wN/Xl1mZBdhUoxheL/m96MLf55asCq4MXIwgVoAXz8BgxfZvL9cqiyCun/fHwW2UaDFCSsnI48kIzuwVPj+/+e7f5/eMvyFOY1RjEFCC8IECCCE9f+/PWH49e7v241/PzEwsoky6nT8E9cGxThDDPOH8/+vVv26dZSBWeD/3/+MH54xCIuhhP/v9/9/vECyW4xZIpVdzYNFVISB4fOf16d/XOv4++IhNGz+v2H8/vE/jyCQCRBALLBy9R/Dj7egtMbKzqhZ9U/aBCgES36MfwSMGK23Mb858u9qDcNvZoZPjxgYDKCZ7dNNhif7/z/f+Pf7a6jpjGxMQsGculEsYpCihoNNxpWJl/3LzqT/3yFG/vh/KJ5JNoBBxAAggKDFDrCOZ3p95P8P1v/S5v+Z8DQA/zI/3/7vh+g/RXNQgv//keFk6N+Xr1BLZyVW+UZuazsGJqSk9/fzl8PWfx5+Qk0crAABxALzItNfUTuCaeY/sHiU9CFUjbAysqEVscAoZGPiAEbEJzS1AAFEWTeDkYmBkRUY1/CYAjnx77d/39/8R8mc/4Cl2d9PD1HCEaiFiQkggFgosf0vAy+j2UbmH+8Ynu/5d3fi/28/waIv/n7c8euxKbuiMLT4+vfr98vtf1/8gaVNBkbRCEYVfwY+aYAAYqRKBxuYBllvFPy6fQDYTITkbCYeO1alZFZpaQaG97+fbf91bfb/37CynFGK0XvLf2ZQwAMEEAsDVQDjv/9/fiDxf/37uvfn1f2/rjOC66r/DChx8YMBWOQBGzUMDAABRJ0uJsvrrX+en4F5HZrYGf7+BRa9IE+jFtIM/98znJ4MUsDAABBAVLCeifEPw6u9/3/8JVrH//9vNjP+BIUWQABRwfp//1n+q5Qwi4qBEhqwWOZzZJEJYuLlQc4gjFLpzKaTGDlYwImRncl81n8OUK0NEECM1BrbYfn1/P/TQ/95DP6LqIFKp6fzfp+bBi1kmDkZXQ/+Z2cDqbkx779s3D9BWYgugACiUtIDVrRskgyK4bDSiYVZQBfYggDnNGBDWeA/GxtEDaNe5X+k4hAggKhmPXqMcOkwKUYyfgQ2Olj/y4bBMjwDst1AABBAjLQbWGP8D02M/xmZcakBCCBa+R6/rXAAEEADPLQEEGAAHnXiBWopm1oAAAAASUVORK5CYII=',
    'settings-2.png': 'iVBORw0KGgoAAAANSUhEUgAAACoAAAAqCAIAAABKoV4MAAAACXBIWXMAAAsTAAALEwEAmpwYAAAABGdBTUEAALGOfPtRkwAAACBjSFJNAAB6JQAAgIMAAPn/AACA6QAAdTAAAOpgAAA6mAAAF2+SX8VGAAAHY0lEQVR42mJ89+4dw8ABgABiYhhQABBAA2w9QAANsPUAAUTA+u9vnr1694N21gMEED7rmf9+vXHz5tWr17/+oZX1AAGEz/r39288evH2zYsHtx99oJH1AAGE03qm/z9vP3r9l0tIiPv/h+cPv9EmAAACCKf1nx/fevT+h5CijraS2K93j+4/+0gL6wECiAmH13/fvvf0H6eUirKcgpqaJB/zxxdPv/35T3XrAQIIu/Wfnt2+/+a7gLyGnCA7h6CssozQt9f3Hj1//5/aDgAIIBYI9R8I/v378/fP71+/vn39/PDu078c4qpqktwsTIyMnGKKCmJP39++coGDUVOIn5eNlZWFlYWZCSjFSKH1AAHE8ubZ4w9fvn779uPZixeMTCwsrGDIxq+gpq4kwsoMMp+RVVBWRekD44v3zx7cefjnNxD++fdPXFyKlZ2dh5dfVFSUi52ZPOsBAojl0c3Lj9/9+M8pICarKivCz83NycXJxcXFwcnJycYE8xwTj5SWsbDyjx/fvgPhN2D4vH3x9OnjVx9/icup8goIkW09QACxGNg48F2/dvvFVw5eUWkFKSEuVqwBysjMxsEFRHwCDP///Pj4+v+Xt69YpBUUVNRUhXjYyA58gABiYmLlUtbSVhHjfH7j7JkbT99/x5++///9+en1/WuXr9/7wyOjoq4qxMtGSdwDBBAo6TGycKrq6jMwXLhx7fSJf4yWOjICHMyM2O3+/OrO5Qs3Hv7jldfSVBPgYsVj9Isbl1/9+PsP5EcudS01ThYsagACCCbGzK6qZ/j/14mrV48f+W/rbCLNhWn/v59v7148f+X+Hz5FfR0NPk588c3M8Ovh46ePP3wGhyWzgIyCghCWcAIIIKR8z8Smpq/H8/fH+6evvmKPgK/vX3348o1dRUMFq93f3jy+eOHW19/gtPL7y6tPX2DG/P/wBVptfnzx4M69J99+Q2UAAgglRP7/+/OTkZlFgJ8Xe/Lj4uZlY2f79O/fP0zJj4+vnb768P3nr6/ePldQVfr65M6vf3BP/Htx68LtPzLf3r1+/fbd+09fn7zWNDbQ5GVnAgggFOt/f33/h5GZj5+fDWb9/z8//zCysoBKGKD17Nw8HGxsvz9/+8kgxIWs8eeHx2dAdn8Bsj+8e3Xh5Cv0GuT9ywtnX8K5rx9dv8bDbaipCBBALKiB+5qBmYVfgIcRlMy+v31w7dKlmy/+CStq6mgoivNzMHNycwOLvHtP36hKCzEzImKInU9Civfmu8/Ep3i+f0wc/zj5AAII1fpPv5iYOPk4/7x/eOPqpRsP33379Y+JieHpzdOvHt6WVtLQUOZhZ2Nj/fPt+59//5mRY5+JVc3I5PnL3W//EmM5m6CskrKKghgnE0AAIaxn/P/ny5c/f35+f3h6161Pn9m5hNSMjJSlRZgY/3989eTOvcc3Tx98xM365+v3vwzff/3+x87MhFqAMPxDTrDMvDLa+uoywrxA175/evva7afvv/6COu7vfwZGRiZmYKEKEEAsyJn61dsPf/8x/f/DJ6NhoiYvwQrN/Iz8YrLGotLvXz66e//x47fv/zJ++/n7Dy8HNCP9+v7t8+ePL549+gRLkYxMHBKGthbqorwcrCyM//8KCgrysp87dfXhmy+/QE78++Xl/etnv70T4QcIIGTr/3KJSmlJyavISbKzYlTEjEyCEgrGYjKKL548fvb531+oVb8+v7l89dbTly///Pn7F1Y+M3JIq6tJCnBDw4eFjYtPXFFe8tHbz1+B1QtQ5M/3j2+eff34mg0ggJACn53PxNYKf6QBq0RhKQVhKaTM8vPr4ydPf/9HVcUlKMSFUm4ysnDz8XKyszEx/PgLyeJ/f//6+/s3QADRop3///9vYNr4j1Ze//n9998/9OIMIIAotZ5LUMrW1sbYUE+QmxVeeP3/9uz+s29//yNZ/vXVizdfvv74DwsMdi5+YQlpOYAAYqHQekZmVmFxSSExCVEBrqMHT3z+Byk83z+4dIHnDzDf8HIw//v95e2zB7fuv/wMq02ZuISkFZWVZIW5AQKIhSrBDWx18fAJcDAyfIb59svzm5d/f34txMMOtP7ru1cv33/++fsfLG2wcgsA/S4jyg4QQCzUivBPb16iVFT//3x9/ejua6xp49/3D2/evH4rzicBEEDUSXr/f3+9c+POt39EKv/388PLZ0+fvvr8GyCAmKgU+Mw8/HyQdikjExs3vwAfBwtKxmNm5+bl4+ZgA4sCW27sStIioszfAAKISoHPwqGqrfOfhe3d138CYpLK8hLfX907ePzSL1h0iMmpGWsr/Pvx5emjh59+MorIKChKCADFAQKIanHPxMGnoasPrBC4uTmAXG4peV6my+/+QjomTIKSMtycHAycHOp8/L//MbKxQu0FCCAWapY3wBQNy/1/GdnlZKS+Pn/9D9xREOXjgGdUNqSqEiCAGGk3tPTn+9dvv/9CShk+Xi6sagACiIWBZoCFk5uPk4AagAAa4LEdgAAaYOsBAmiArQcIMAAM7+OMMd1d/AAAAABJRU5ErkJggg=='
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* MonomerManager class
* @class org.helm.webeditor.MonomerManager
* Recommended database schema:
* <pre>
* **********************************************************
* create table HELMMonomers
* (
* id bigint not null identity(1, 1) primary key,
* Symbol varchar(256) not null,
* Name varchar(256) not null,
* NaturalAnalog varchar(256),
* SMILES varchar(max),
* PolymerType varchar(256) not null,
* MonomerType varchar(256),
* Molfile varchar(max),
* Hashcode varchar(128),
* R1 varchar(256),
* R2 varchar(256),
* R3 varchar(256),
* R4 varchar(256),
* R5 varchar(256),
* Author nvarchar(256),
* CreatedDate DateTime default getdate()
* );
* **********************************************************
* </pre>
* JSON Schema
* <pre>
* {
*     id: 69,                     // monomer internal ID
*     symbol: 'Alexa',            // monomer symbol
*     name: 'Alexa Fluor 488',    // monomer long name
*     naturalanalog: null,        // natural analog
*     smiles: null,               // smiles
*     polymertype: 'CHEM',        // polymer type: CHEM, SUGAR, LINKER, BASE, AA
*     monomertype: null,          // momer type: Backbone, Branch, null
*     molfile: null,              // molfile of monomer, plain text, not BASE64 encoded or compressed
*     r1: 'X',                    // cap for R1 
*     r2: null,                   // cap for R2
*     r3: null,                   // cap for R3
*     r4: null,                   // cap for R4
*     r5: null,                   // cap for R5
*     author: null,               // monomer author
*     createddate: null           // monomer created date
* }
* </pre>
**/
org.helm.webeditor.MonomerManager = scil.extend(scil._base, {
    /**
    * @constructor MonomerManager
    * @param {DOM} parent - The parent element to host Monomer Manager
    * @bio {dict} options - options on how to render the App
    * <pre>
    * ajaxurl: {string} The service url for the ajax
    * <b>Example:</b>
    *     <div id="div1" style="margin: 5px; margin-top: 15px"></div>
    *     <script type="text/javascript">
    *       scil.ready(function () {
    *         new org.helm.webeditor.MonomerManager("div1", { ajaxurl: "../service/ajaxtool/post?cmd=" });
    *       });
    *     </script>
    * </pre>
    **/
    constructor: function (parent, options) {
        if (typeof (parent) == "string")
            parent = scil.byId(parent);
        this.options = options == null ? {} : options;
        scil.Page.ajaxurl = this.options.ajaxurl;
        this.init(parent);
    },

    /**
    * Initialize the manager (internal use)
    * @function init
    */
    init: function (parent) {
        var me = this;

        this.page = new scil.Page(parent);

        var me = this;
        var ispro = scil.helm != null && JSDraw2.Security.kEdition != "Lite";
        this.buttons = ["-"];
        this.buttons.push({ type: "input", key: "symbol", labelstyle: { fontSize: "90%" }, label: "Symbol", styles: { width: 100 }, autosuggesturl: this.options.ajaxurl + "helm.monomer.suggest", onenter: function () { me.refresh(); }, onchange: function () { me.clearFilterValue("name"); } });
        if (!ispro)
            this.buttons.push({ type: "input", key: "name", labelstyle: { fontSize: "90%" }, label: "Name", styles: { width: 100 }, autosuggesturl: this.options.ajaxurl + "helm.monomer.suggest", onenter: function () { me.refresh(); }, onchange: function () { me.clearFilterValue("symbol"); } });
        this.buttons.push({ type: "select", key: "polymertype", labelstyle: { fontSize: "90%" }, items: org.helm.webeditor.MonomerManager.getPolymerTypes(), label: "Polymer Type", styles: { width: 100 }, onchange: function () { me.refresh(); } });
        this.buttons.push({ type: "select", key: "monomertype", labelstyle: { fontSize: "90%" }, items: org.helm.webeditor.MonomerManager.getMonomerTypes(), label: "Monomer Type", styles: { width: 100 }, onchange: function () { me.refresh(); } });
        this.buttons.push({ type: "select", key: "countperpage", labelstyle: { fontSize: "90%" }, label: "Count", items: ["", 10, 25, 50, 100], onchange: function () { me.refresh(); } });

        if (typeof (JSDrawServices) != "undefined" && JSDrawServices.url != null) {
            this.buttons.splice(0, 0, "-");
            this.buttons.splice(0, 0, { type: "a", src: scil.Utils.imgSrc("img/open.gif"), title: "Import Monomers", onclick: function () { me.uploadFile(true); } });
            this.buttons.splice(0, 0, "-");
            this.buttons.splice(0, 0, { type: "a", src: scil.Utils.imgSrc("img/save.gif"), title: "Export Monomers", items: ["JSON", "SDF"], onclick: function (cmd) { me.exportFile(cmd); } });
        }

        if (org.helm.webeditor.MonomerManager.onGetButtons != null)
            org.helm.webeditor.MonomerManager.onGetButtons(this.buttons, this);

        var columns = org.helm.webeditor.MonomerManager.columns;
        var fields = org.helm.webeditor.MonomerManager.fields;
        fields.polymertype.items = org.helm.webeditor.MonomerManager.getPolymerTypes();
        fields.monomertype.items = org.helm.webeditor.MonomerManager.getMonomerTypes();
        for (var i = 1; i <= 5; ++i)
            fields["r" + i].items = org.helm.webeditor.MonomerManager.caps;

        if (ispro) {
            scil.helm.MonomerManager.ajaxurl = this.options.ajaxurl;
            columns.versions = { label: "Versions", type: "html", render: function (v, values) { return v > 0 ? "<a href='javascript:scil.helm.MonomerManager.showVersions(" + values.id + ")'>" + v + "</a>" : null; } };
        }
        else {
            columns.aliases = null;
            fields.aliases = null;
        }

        this.monomers = this.page.addForm({
            caption: "Monomer List",
            key: "id",
            object: "helm.monomer",
            imagewidth: 30,
            buttons: this.buttons,
            onbeforerefresh: function (args) { me.onbeforerefresh(args); },
            onbeforesave: function (data, args, form) { return me.onbeforesave(data, args, form); },
            onshowform: function (dlg, args, action) { dlg.form.fields.symbol.disabled = action != "create"; },
            columns: columns,
            savedoc: true,
            formcaption: "Monomer",
            fields: fields
        });

        this.page.addForm({
            caption: "Monomer",
            type: "form",
            object: "helm.monomer",
            fields: org.helm.webeditor.MonomerManager.fields
        }, this.monomers);

        this.monomers.refresh();
    },

    clearFilterValue: function (key) {
        for (var i = 0; i < this.buttons.length; ++i) {
            if (this.buttons[i].key == key) {
                this.buttons[i].b.value = "";
                break;
            }
        }
    },

    onbeforesave: function (data, args, form) {
        if (data.polymertype != "CHEM" && scil.Utils.isNullOrEmpty(data.naturalanalog)) {
            scil.Utils.alert("Natural Analog cannot be blank");
            return false;
        }

        data.molfile = form.fields.molfile.jsd.getMolfile();

        // check R caps
        var ratoms = {};
        var atoms = form.fields.molfile.jsd.m.atoms;
        for (var i = 0; i < atoms.length; ++i) {
            var a = atoms[i];
            if (a.elem == "R") {
                var r = (a.alias == null ? "R" : a.alias);
                if (ratoms[r.toLowerCase()] != null) {
                    scil.Utils.alert("The R cannot be used twice: " + r);
                    return false;
                }
                ratoms[r.toLowerCase()] = r;
            }
        }

        for (var r in ratoms) {
            var cap = data[r];
            if (scil.Utils.isNullOrEmpty(cap)) {
                scil.Utils.alert("The cap of " + ratoms[r] + " is not defined yet");
                return false;
            }
        }

        for (var i = 1; i <= 5; ++i) {
            var r = "r" + i;
            if (!scil.Utils.isNullOrEmpty(data[r.toLowerCase()]) && ratoms[r.toLowerCase()] == null) {
                scil.Utils.alert("R" + i + " is defined, but not drawn in the structure");
                return false;
            }
        }
    },

    /**
    * Refresh the list (internal use)
    * @function refresh
    */
    refresh: function (view) {
        this.monomers.refresh();
    },

    /**
    * Event handler before refreshing (internal use)
    * @function onbeforerefresh
    */
    onbeforerefresh: function (args) {
        scil.Form.getButtonValuesByKey(this.buttons, org.helm.webeditor.MonomerManager.kFilters, args);
    },

    /**
    * Import from file (internal use)
    * @function uploadFile
    */
    uploadFile: function (duplicatecheck) {
        scil.Utils.uploadFile("Import Monomer Library", "Select HELM monomer xml file, json or SDF file (" + (duplicatecheck ? "with" : "without") + " duplicate check)", this.options.ajaxurl + "helm.monomer.uploadlib",
            function (ret) { scil.Utils.alert2(ret.n + " monomers are imported.\r\n" + (ret.errors == null ? "" : ret.errors)); }, { duplicatecheck: duplicatecheck });
    },

    exportFile: function (ext) {
        window.open(this.options.ajaxurl.replace("/post?", "/get?") + "helm.monomer.savefile&wrapper=raw&ext=" + ext, "_blank");
    }
});

scil.apply(org.helm.webeditor.MonomerManager, {
    caps: ["", "H", "OH", "X"],
    capsmiles: { "H": "[H]", "OH": "O", "X": "[X]" },
    kFilters: ["status", "polymertype", "monomertype", "status", "symbol", "name", "countperpage"],

    showVersions: function (id) {
        var me = this;
        scil.Utils.ajax(this.ajaxurl + "helm.monomer.versions", function (ret) {
            me.showVersions2(ret);
        }, { id: id });
    },

    showVersions2: function (ret) {
        if (this.versionDlg == null) {
            var fields = {
                table: {
                    type: "table", viewonly: true, columns: {
                        versionid: { label: "Version ID" },
                        dt: { label: "Date", type: "date" },
                        user: { label: "User" },
                        action: { label: "Action" }
                    }
                }
            };
            this.versionDlg = scil.Form.createDlgForm("Versions", fields, null, { hidelabel: true });
        }

        var cur = ret[ret.length - 1];
        cur.versionid = "Current";
        ret.splice(ret.length - 1, 1);
        ret.splice(0, 0, cur);

        this.versionDlg.show();
        this.versionDlg.form.setData({ table: ret });
        this.versionDlg.moveCenter();
    },

    /**
    * Tool function (internal use)
    * @function getValueByKey
    */
    getValueByKey: function (list, key) {
        for (var i = 0; i < list.length; ++i) {
            if (list[i].key == key)
                return list[i].b.value;
        }
        return null;
    },

    /**
    * List polymer types (internal use)
    * @function getPolymerTypes
    */
    getPolymerTypes: function () {
        return ["", "RNA", "CHEM", "PEPTIDE"];
    },

    /**
    * List monomer types (internal use)
    * @function getMonomerTypes
    */
    getMonomerTypes: function () {
        return ["", "Backbone", "Branch", "Undefined"]
    },

    /**
    * List of statuses (internal use)
    * @function getStatuses
    */
    getStatuses: function () {
        return ["", "New", "Approved", "Retired"]
    },

    columns: {
        id: { type: "hidden", iskey: true },
        symbol: { label: "Symbol", width: 100 },
        aliases: { label: "Aliases", width: 100 },
        name: { label: "Name", width: 200 },
        naturalanalog: { label: "Natural Analog", width: 100 },
        polymertype: { label: "Polymer Type", width: 100 },
        monomertype: { label: "Monomer Type", width: 100 },
        r1: { label: "R1", width: 50 },
        r2: { label: "R2", width: 50 },
        r3: { label: "R3", width: 50 },
        author: { label: "Author", width: 100 },
        versions: null,
        createddate: { label: "Created Date", type: "date", width: 100 }
    },

    fields: {
        id: { type: "hidden" },
        symbol: { label: "Symbol", required: true },
        aliases: { label: "Aliases", width: 800 },
        name: { label: "Name", required: true, width: 800 },
        polymertype: { label: "Polymer Type", required: true, type: "select", items: null, width: 100 },
        monomertype: { label: "Monomer Type", required: true, type: "select", items: null, width: 100 },
        naturalanalog: { label: "Natural Analog", required: true, width: 100 },
        author: { label: "Author", width: 100 },
        smiles: { label: "SMILES", width: 800, viewonly: true },
        molfile: { label: "Structure", type: "jsdraw", width: 800, height: 300 },
        r1: { label: "R1", type: "select", items: null },
        r2: { label: "R2", type: "select", items: null },
        r3: { label: "R3", type: "select", items: null },
        r4: { label: "R4", type: "select", items: null },
        r5: { label: "R5", type: "select", items: null }
    }
});


org.helm.webeditor.MonomerLibApp = org.helm.webeditor.MonomerManager;
scil.helm.MonomerLibApp = org.helm.webeditor.MonomerManager;﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* RuleSet class
* @class org.helm.webeditor.RuleSet
*/
org.helm.webeditor.RuleSet = {
    kApplyAll: false,

    rules: [
        { id: 1, category: "Demo", name: "Replace base A with U", description: "", script: "function(plugin) {var n = plugin.replaceMonomer(org.helm.webeditor.HELM.BASE, 'A', 'U');return n > 0;}" },
        { id: 2, category: "Demo", name: "Replace base A with G", description: "", script: "function(plugin) {var n = plugin.replaceMonomer(org.helm.webeditor.HELM.BASE, 'A', 'G');return n > 0;}" },
        { id: 3, category: "Test", name: "Replace base A with T", description: "", script: "function(plugin) {var n = plugin.replaceMonomer(org.helm.webeditor.HELM.BASE, 'A', 'T');return n > 0;}" }
    ],

    loadDB: function (list) {
        this.rules = list;
    },

    favorites: new scil.Favorite("ruleset"),

    saveTextDB: function (url) {
        var cols = ["id", "name", "description", "script", "author", "category"];

        var n = 0;
        var ret = "";
        for (var i = 0; i < this.rules.length; ++i) {
            var r = this.rules[i];
            var s = "";
            for (var k = 0; k < cols.length; ++k)
                s += (k > 0 ? "|" : "") + r[cols[k]];
            ret += JSDraw2.Base64.encode(s) + "\n";
            ++n;
        }

        ret = n + "\n" + ret;
        if (url == null)
            return ret;

        var args = { client: "jsdraw", wrapper: "none", filename: "rules.txt", directsave: 1, contents: ret };
        scil.Utils.post(url, args, "_blank");
    },

    addFavorite: function (e) {
        var img = e.srcElement || e.target;
        var tr = scil.Utils.getParent(img, "TR");
        var id = tr.getAttribute("ruleid");

        var f = img.getAttribute("star") != "1";
        if (f) {
            img.setAttribute("star", "1");
            img.src = scil.Utils.imgSrc("img/star.png");
        }
        else {
            img.setAttribute("star", "");
            img.src = scil.Utils.imgSrc("img/star0.png");
        }

        this.favorites.add(id, f);
    },

    filterRules: function (tbody, s, category) {
        s = scil.Utils.trim(s).toLowerCase();
        var list = tbody.childNodes;
        for (var i = 0; i < this.rules.length; ++i) {
            var r = this.rules[i];
            var tr = list[i + 1];
            if ((s == "" || r.name.toLowerCase().indexOf(s) >= 0) && (scil.Utils.isNullOrEmpty(category) || category == r.category))
                tr.style.display = "";
            else
                tr.style.display = "none";
        }
    },

    listRules: function (mex, apply, applyall) {
        var div = mex.divRule;
        scil.Utils.removeAll(div);

        var me = this;
        var tbody = scil.Utils.createTable(div, 0, 0, { width: "100%" });
        var tr = scil.Utils.createElement(tbody, "tr", null, { background: "#eee", display: (this.kApplyAll ? "" : "none") });
        var chk = scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "checkbox");
        scil.Utils.createButton(scil.Utils.createElement(tr, "td", null, { textAlign: "right", padding: "3px 3px 3px 0" }, { colSpan: 3 }), this.createApplyAll("Apply All", applyall, tbody));
        scil.connect(chk, "onclick", function () { me.checkAll(tbody); });

        var k = 1;
        var list = [];
        for (var i = 0; i < this.rules.length; ++i) {
            var r = this.rules[i];
            var fav = this.favorites.contains(r.id);
            if (this.favorites.contains(r.id))
                this.listOneRule(mex, tbody, r, ++k, apply, true);
            else
                list.push(r);
        }

        for (var i = 0; i < list.length; ++i)
            this.listOneRule(mex, tbody, list[i], ++k, apply);

        return tbody;
    },

    listOneRule: function (mex, tbody, r, i, apply, fav) {
        var me = this;
        var tr = scil.Utils.createElement(tbody, "tr", null, { background: i % 2 == 1 ? "#eee" : null }, { ruleid: r.id });
        scil.Utils.createElement(scil.Utils.createElement(tr, "td"), "checkbox", null, { display: (this.kApplyAll ? "" : "none"), width: "1%" });

        var td = scil.Utils.createElement(tr, "td");
        scil.Utils.createElement(td, "img", null, { /*width: "1%"*/
        }, { star: (fav ? 1 : null), src: scil.Utils.imgSrc("img/star" + (fav ? "" : "0") + ".png") }, function (e) { me.addFavorite(e); mex.listRules(); });

        td = scil.Utils.createElement(tr, "td", null, { width: "99%" });
        this.listOneRule2(td, r, apply, i);
    },

    listOneRule2: function (td, rule, fun, i) {
        var s = rule.name;
        if (scil.Utils.isNullOrEmpty(s))
            s = rule.description;
        if (s.length > 50)
            s = s.substr(0, 47) + "...";

        var tbody = scil.Utils.createTable(td, 0, 0, { width: "100%" });
        var tr = scil.Utils.createElement(tbody, "tr");
        scil.Utils.createElement(tr, "td", "[" + rule.id + "] " + s, { padding: "3px 0 3px 0" }, { title: rule.description });
        var button = scil.Utils.createElement(scil.Utils.createElement(tr, "td", null, { textAlign: "right" }), "button", JSDraw2.Language.res("Apply"), { display: "none" });

        var me = this;
        scil.connect(button, "onclick", function () { fun(rule.script); });
        scil.connect(td, "onmouseover", function (e) { button.style.display = ""; });
        scil.connect(td, "onmouseout", function (e) { button.style.display = "none"; });
    },

    checkAll: function (tbody) {
        var nodes = tbody.childNodes;
        var f = nodes[0].childNodes[0].childNodes[0].checked;
        for (var i = 1; i < nodes.length; ++i) {
            var tr = nodes[i];
            tr.childNodes[0].childNodes[0].checked = f;
        }
    },

    createApplyAll: function (label, fun, tbody) {
        return {
            label: label, type: "a", onclick: function (e) {
                var list = [];
                var nodes = tbody.childNodes;
                for (var i = 1; i < nodes.length; ++i) {
                    var tr = nodes[i];
                    if (tr.childNodes[0].childNodes[0].checked)
                        list.push(parseInt(tr.getAttribute("ruleid")));
                }

                if (list.length == 0)
                    scil.Utils.alert("No rule selected");
                else
                    fun(list);
            }
        };
    },

    applyRules: function (plugin, ruleids) {
        if (ruleids.length == 0)
            return;

        var list = [];
        for (var i = 0; i < ruleids.length; ++i) {
            for (var k = 0; k < this.rules.length; ++k) {
                var r = this.rules[k];
                if (ruleids[i] == r.id) {
                    list.push(r);
                    break;
                }
            }
        }

        var args = { plugin: plugin, n: list.length, changed: 0, list: list, cloned: plugin.jsd.clone() };
        this._applyNextRule(args);
    },

    applyRule: function (plugin, script) {
        var list = [{ script: script, name: null}];
        var args = { plugin: plugin, n: list.length, changed: 0, list: list, cloned: plugin.jsd.clone() };
        this._applyNextRule(args);
    },

    _applyNextRule: function (args) {
        if (args.list.length == 0)
            return;

        var me = this;

        // get the first rule 
        var rule = args.list[0];
        args.list.splice(0, 1);

        // callback function when the rule is applied
        var callback = function (f, error) {
            if (error != null) {
                // some rule failed
                scil.Utils.alert(error);
                args.plugin.jsd.restoreClone(args.cloned);
                return;
            }

            if (f)
                ++args.changed; // structure changed

            if (args.list.length > 0) {
                // continue to apply the next rule
                me._applyNextRule(args);
                return;
            }

            // all rules are applied
            if (args.changed > 0) {
                args.plugin.jsd.pushundo(args.cloned);
                args.plugin.jsd.refresh(true);
                scil.Utils.alert((args.n > 1 ? "Rules" : "Rule") + " applied successfully!");
            }
            else {
                scil.Utils.alert((args.n > 1 ? "Rules" : "Rule") + " applied, but nothing changed!");
            }
        };
        this._applyOneRule(args.plugin, rule.script, rule.name, callback);
    },

    _applyOneRule: function (plugin, script, name, callback) {
        var rulefun = null;
        if (typeof (script) == "string")
            rulefun = scil.Utils.eval(script);
        else if (typeof (script) == "function")
            rulefun = script;

        var f = false;
        var error = null;
        if (rulefun == null) {
            error = "Error: Invalid rule function: " + name;
        }
        else {
            try {
                f = rulefun(plugin);
            }
            catch (e) {
                error = "Error: " + (name == null ? "" : name) + "\n---------------\n" + e.message + "\n---------------\n" + e.stack;
            }
        }

        callback(f, error);
    }
};﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

/**
* RuleManager class
* @class org.helm.webeditor.RuleManager
* Recommended database schema:
* <pre>
* **********************************************************
* create table HELMRules
* (
* id bigint not null identity(1, 1) primary key,
* Category nvarchar(256),
* Name nvarchar(512) not null,
* Script nvarchar(max),
* Description nvarchar(max),
* Author nvarchar(256),
* CreatedDate DateTime default getdate()
* );
* **********************************************************
* </pre>
* JSON Schema
* <pre>
* {
*     id: 3,                                          // rule internal ID
*     name: 'Replace base A with U',                  // rule long name
*     script: '	\nfunction(plugin) {\n // ... \n}',   // rule script
*     description: null,                              // rule full description
*     author: null,                                   // rule author
*     createddate: null,                              // rule created date
*     category: null                                  // rule category
* }
* </pre>
**/
org.helm.webeditor.RuleManager = scil.extend(scil._base, {
    /**
    * @constructor RuleManager
    * @param {DOM} parent - The parent element to host the Ruleset Manager
    * @bio {dict} options - options on how to render the App
    * <pre>
    * ajaxurl: {string} The service url for the ajax
    * <b>Example:</b>
    *     <div id="div1" style="margin: 5px; margin-top: 15px"></div>
    *     <script type="text/javascript">
    *       scil.ready(function () {
    *         new org.helm.webeditor.RuleManager("div1", { ajaxurl: "../service/ajaxtool/post?cmd=" });
    *       });
    *     </script>
    * </pre>
    **/
    constructor: function (parent, options) {
        if (typeof (parent) == "string")
            parent = scil.byId(parent);
        this.options = options == null ? {} : options;
        scil.Page.ajaxurl = this.options.ajaxurl;
        this.init(parent);
    },

    init: function (parent) {
        var me = this;

        this.page = new scil.Page(parent);

        var me = this;
        this.buttons = [
            "-",
            { type: "select", key: "category", labelstyle: { fontSize: "90%" }, items: org.helm.webeditor.RuleManager.categories, label: "Category", styles: { width: 100 }, onchange: function () { me.refresh(); } },
            { type: "select", key: "countperpage", labelstyle: { fontSize: "90%" }, label: "Count", items: ["", 10, 25, 50, 100], onchange: function () { me.refresh(); } }
        ];

        var fields = org.helm.webeditor.RuleManager.getFields();
        fields.script.button = [{ label: "Test Script", onclick2: function (field) { me.testscript(field); } },
            { label: "Test Applying", onclick2: function (field, form) { me.testapplying(field, form); } }
        ];
        fields.test = { label: "Test Structure", type: "jsdraw", width: 800, height: 300 };
        this.rules = this.page.addForm({
            caption: "Rule Set",
            key: "id",
            object: "helm.rule",
            buttons: this.buttons,
            onbeforerefresh: function (args) { me.onbeforerefresh(args); },
            onbeforesave: function (data, args, form) { data.test = null; },
            columns: {
                id: { label: "ID", width: 50, iskey: true },
                category: { label: "Category", width: 60 },
                name: { label: "Name", width: 100 },
                description: { label: "Description", width: 200 },
                author: { label: "Author", width: 100 },
                createddate: { label: "Created Date", type: "date", width: 100 }
            },
            formcaption: "Rule",
            fields: fields,
            defaultvalues: { script: "function(plugin) {\n\n}\n\n//function(plugin) { \n//    scil.Utils.ajax('http://SERVER/youerservice', function(ret) {\n//        plugin.setHelm(ret.new_helm);\n//    });\n//} " }
        });

        this.page.addForm({
            caption: "Rule",
            type: "form",
            object: "helm.rule",
            fields: org.helm.webeditor.RuleManager.getFields()
        }, this.rules);

        this.rules.refresh();
    },

    refresh: function (view) {
        this.rules.refresh();
    },

    onbeforerefresh: function (args) {
        scil.Form.getButtonValuesByKey(this.buttons, ["category", "countperpage"], args);
    },

    testapplying: function (field, form) {
        var plugin = form.fields.test.jsd.helm;
        plugin.applyRule(field.value);
    },

    testscript: function (field) {
        var s = field.value;
        if (scil.Utils.trim(s) == "")
            return;

        try {
            eval("var __fun=" + s);
            if (typeof (__fun) == "function")
                scil.Utils.alert("Looks good!");
            else
                scil.Utils.alert("It should be a Javascript function, like this: \nfunction(plugin) {\n //... \n}");
        }
        catch (e) {
            scil.Utils.alert(e.message);
        }
    }
});

scil.apply(org.helm.webeditor.RuleManager, {
    categories: ["", "General"],

    getFields: function () {
        return {
            id: { label: "ID", viewonly: true },
            category: { label: "Category", width: 200, type: "select", items: this.categories },
            name: { label: "Name", width: 800 },
            description: { label: "Description", type: "textarea", width: 800, height: 40 },
            author: { label: "Author", width: 100 },
            script: { label: "Javascript", type: "textarea", width: 800, height: 160 }
        }
    }
});


org.helm.webeditor.RuleSetApp = org.helm.webeditor.RuleManager;﻿/*******************************************************************************
* Copyright (C) 2018, The Pistoia Alliance
* Created by Scilligence, built on JSDraw.Lite
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the 
* following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

org.helm.webeditor.Adapter = {
    init: function (url) {
        org.helm.webeditor.Adapter.url = url; // http://localhost:8080/HELMMonomerService-master/webService/service/library

        scil.Utils.onAjaxCallback = org.helm.webeditor.Adapter.onAjaxCallback;
        scil.Utils.onajaxcall = org.helm.webeditor.Adapter.onajaxcall;
    },

    startApp: function (div, options) {
        org.helm.webeditor.Adapter.init(options.url);

        if (options.onValidateHelm == null)
            options.onValidateHelm = this.onValidateHelm;
        if (options.onCleanUpStructure == null)
            options.onCleanUpStructure = this.onCleanUpStructure;
        if (options.onMonomerSmiles == null)
            org.helm.webeditor.Monomers.onMonomerSmiles = this.onMonomerSmiles;
        else
            org.helm.webeditor.Monomers.onMonomerSmiles = options.onMonomerSmiles;

        scil.Utils.ajax(org.helm.webeditor.Adapter.url + "/monomer/ALL", function (ret) {
            org.helm.webeditor.Adapter.startApp2(div, options);
        });

        scil.Utils.ajax(org.helm.webeditor.Adapter.url + "/rule", function (ret) {
            org.helm.webeditor.RuleSet.rules = ret;
        });
    },

    onAjaxCallback: function (ret) {
        if (scil.Utils.isNullOrEmpty(ret))
            return;

        // rules:
        // TODO ...

        // monomers:
        var list = ret.monomers;
        if (list == null) {
            org.helm.webeditor.Adapter.toHWE(ret, ret);
            return;
        }

        if (ret.limit == null || ret.limit == 0) {
            org.helm.webeditor.Adapter.loadMonomers(ret);
        }
        else {
            var rows = [];
            for (var i = 0; i < list.length; ++i) {
                var x = list[i];
                var m = { symbol: x.symbol, name: x.name, createddate: x.createDate };
                org.helm.webeditor.Adapter.toHWE(m, x);
                rows.push(m);
            }
            if (ret.offset != null && ret.limit > 0) {
                ret.page = ret.offset / ret.limit + 1;
                var mod = ret.total % ret.limit;
                ret.pages = (ret.total - mod) / ret.limit + (mod > 0 ? 1 : 0);
            }
            ret.rows = rows;
        }
    },

    onajaxcall: function (args, opts) {
        args.headers = { 'Accept': 'application/json', 'Content-Type': 'application/json;charset=utf-8' };
        switch (args.url) {
            case "helm.monomer.load":
                args.url = org.helm.webeditor.Adapter.url + "/monomer/" + args.content.id;
                opts.verb = "get";
                break;
            case "helm.monomer.list":
                var limit = args.content.countperpage;
                if (!(limit > 0))
                    limit = 10;
                var page = args.content.page;
                if (!(page > 0))
                    page = 1;
                var offset = (page - 1) * limit;
                var pt = scil.Utils.isNullOrEmpty(args.content.polymertype) ? "ALL" : args.content.polymertype;
                args.url = org.helm.webeditor.Adapter.url + "/monomer/" + pt + "?limit=" + limit + "&offset=" + offset;
                if (!scil.Utils.isNullOrEmpty(args.content.symbol))
                    args.url += "&filter=" + escape(args.content.symbol) + "&filterField=symbol";
                else if (!scil.Utils.isNullOrEmpty(args.content.name))
                    args.url += "&filter=" + escape(args.content.name);// + "&filterField=name";
                opts.verb = "get";
                break;
            case "helm.monomer.save":
                if (scil.Utils.isNullOrEmpty(args.content.id)) // new monomer
                    args.content.id = args.content.polymertype + "/" + args.content.symbol;
                args.url = org.helm.webeditor.Adapter.url + "/monomer/" + args.content.id;
                opts.verb = "put";
                args.content.id = null;
                org.helm.webeditor.Adapter.fromHWE(args.content);
                args.postData = scil.Utils.json2str(args.content, null, true);
                delete args.content;
                break;
            case "helm.monomer.del":
                args.url = org.helm.webeditor.Adapter.url + "/monomer/" + args.content.id;
                opts.verb = "del";
                args.content = {};
                break;

            case "helm.rule.load":
                args.url = org.helm.webeditor.Adapter.url + "/rule/" + args.content.id;
                opts.verb = "get";
                break;
            case "helm.rule.list":
                var limit = args.content.countperpage;
                if (!(limit > 0))
                    limit = 10;
                var page = args.content.page;
                if (!(page > 0))
                    page = 1;
                var offset = (page - 1) * limit;
                args.url = org.helm.webeditor.Adapter.url + "/rule?limit=" + limit + "&offset=" + offset;
                opts.verb = "get";
                break;
            case "helm.rule.save":
                args.url = org.helm.webeditor.Adapter.url + "/rule";
                opts.verb = "put";
                //args.content.id = null;
                org.helm.webeditor.Adapter.fromHWE(args.content);
                args.postData = scil.Utils.json2str(args.content, null, true);
                delete args.content;
                break;
            case "helm.rule.del":
                args.url = org.helm.webeditor.Adapter.url + "/rule/" + args.content.id;
                opts.verb = "del";
                args.content = {};
                break;

            default:
                if (opts.verb == null)
                    opts.verb = "get";
                break;
        }
        opts.ignoresucceedcheck = true;
    },


    startApp2: function (div, options) {
        org.helm.webeditor.ambiguity = options.ambiguity;
        new scil.helm.AppToolbar(options.toolbarholder, "helm/img/", options.toolbarbuttons);
        app = new scil.helm.App(div, options);
    },

    loadMonomers: function (ret) {
        org.helm.webeditor.Monomers.clear();
        var list = ret.monomers;
        for (var i = 0; i < list.length; ++i) {
            var x = list[i];
            var m = { id: x.symbol, n: x.name, na: x.naturalAnalog, type: x.polymerType, mt: x.monomerType, m: x.molfile };

            m.at = {};
            var rs = 0;
            if (x.rgroups != null) {
                for (var k = 0; k < x.rgroups.length; ++k) {
                    var r = x.rgroups[k];
                    m.at[r.label] = r.capGroupName;
                    ++rs;
                }
            }
            m.rs = rs;

            org.helm.webeditor.Monomers.addOneMonomer(m);
        }
    },

    toHWE: function (m, ret) {
        if (ret.polymerType == null)
            return;

        m.id = ret.polymerType + "/" + ret.symbol;
        m.naturalanalog = ret.naturalAnalog;
        m.polymertype = ret.polymerType;
        m.monomertype = ret.monomerType;
        m.author = ret.author;

        if (ret.rgroups == null)
            return;
        for (var k = 0; k < ret.rgroups.length; ++k) {
            var r = ret.rgroups[k];
            m[scil.helm.symbolCase(r.label)] = r.capGroupName;
        }
    },

    fromHWE: function (ret) {
        ret.naturalAnalog = ret.naturalanalog;
        ret.polymerType = ret.polymertype;
        ret.monomerType = ret.monomertype;

        var rgroups = [];
        for (var i = 1; i < 5; ++i) {
            if (ret["r" + i] != null) {
                var cap = ret["r" + i];
                var smiles = "[*:" + i + "]" + org.helm.webeditor.MonomerManager.capsmiles[cap];
                rgroups.push({ alternateId: "R" + i + "-" + ret["r" + i], label: "R" + i, capGroupName: cap, capGroupSMILES: smiles });
            }
        }
        ret.rgroups = rgroups;
    },

    onValidateHelm: function (me) {
        var url = me.options.validateurl;
        if (scil.Utils.isNullOrEmpty(url)) {
            scil.Utils.alert("The validation url is not configured yet");
            return;
        }

        me.setNotationBackgroundColor("white");
        var helm = scil.Utils.getInnerText(me.notation);
        if (scil.Utils.isNullOrEmpty(helm))
            return;

        scil.Utils.ajax(url,
            function (ret) { me.setNotationBackgroundColor(ret.Validation == "valid" ? "#9fc" : "#fcf"); },
            { HELMNotation: helm },
            { onError: function (data) { me.setNotationBackgroundColor("#fcf"); }, verb: "post", headers: { Accept: "application/json", "Content-Type": "application/x-www-form-urlencoded" }
            });
    },

    onCleanUpStructure: function (mol, me) {
        //        scil.Utils.ajax(me.options.cleanupurl, function (ret) {
        //            me.structureview.setMolfile(ret == null ? null : ret.output);
        //        }, { input: mol.getMolfile(), inputformat: "mol", outputformat: "mol" });
        var url = me.options.cleanupurl;
        var molfile = mol.getMolfile();
        scil.Utils.ajax(url,
            function (ret) { m.m = ret.Molfile; },
            { SMILES: molfile },
            { verb: "post", headers: { Accept: "application/json", "Content-Type": "application/x-www-form-urlencoded" }
            });
    },

    onMonomerSmiles: function (m, smiles) {
        var url = org.helm.webeditor.Monomers.cleanupurl;
        scil.Utils.ajax(url,
            function (ret) { m.m = ret.Molfile; },
            { SMILES: smiles },
            { verb: "post", headers: { Accept: "application/json", "Content-Type": "application/x-www-form-urlencoded" }
            });
    }
};