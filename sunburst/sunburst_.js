const width = document.getElementById("container").offsetWidth * 0.95,
    height = 800
radius = Math.min(width, height) / 2;

// Couleurs de base pour les 9 entreprises sélectionnées
var positionColors = ["#da1d23", "#ebb40f", "#187a2f", "#0aa3b5", "#c94930", "#ad213e", "#a87b64", "#e65832", "#da0d68"];

// Rotation de la teinte pour le niveau 3 (catégorie). HSL = Hue/Saturation/Lightness = Teinte/Saturation/Valeur
var rotationColors = {
    "Académies, Fondation, sociétés savantes, organismes de conseils": 5,
    "Association d'étudiants": 10,
    "Association professionnel de santé": 15,
    "Association usager de santé": 20,
    "Editeur de logiciel": 25,
    "Etablissement de santé": 30,
    "Etudiant": 35,
    "Personnes morales assurant la formation initiale ou continue des professionnels de santé": 40,
    "Presse et média": 45,
    "Professionnel de santé": 50,
    "Vétérinaire": 55,
    "Vétérinaire Personne Morale": 60
};

const vis = d3.select("#chart").append("svg")
    .attr("class", "svg")
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("id", "chart-container")
    .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

// Pour ajouter les patterns
vis.append("defs")
    .attr("id", "defs");

d3.text("d3js/sunburst-chart/transparence_data.csv").then(function (raw) {
    let dsv = d3.dsvFormat(';');
    let data = dsv.parse(raw);
    let json = buildHierarchy(data);
    addTextElement();
    createVisualization(json);
    buildMinChart(json);
});

// Ajout du texte au milieu du cercle
function addTextElement() {
    var textGroup = vis.append("g");

    textGroup.append("text")
        .attr("id", "entreprise")
        .attr("y", -100)
        .attr("class", "entreprise")
        .attr("text-anchor", "middle");

    textGroup.append("text")
        .attr("id", "type-amount")
        .attr("y", -80)
        .attr("class", "type-amount")
        .attr("text-anchor", "middle");
    textGroup.append("text")
        .attr("id", "category-amount")
        .attr("y", -60)
        .attr("class", "category-amount")
        .attr("text-anchor", "middle");
    textGroup.append("text")
        .attr("id", "amount")
        .attr("class", "amount")
        .attr("text-anchor", "middle");
}

function createVisualization(json) {
    var arc = d3.arc()
        .startAngle(function (d) {
            return d.x0;
        })
        .endAngle(function (d) {
            return d.x1;
        })
        .innerRadius(function (d) {
            return Math.sqrt(d.y0);
        })
        .outerRadius(function (d) {
            return Math.sqrt(d.y1);
        });

    var partition = d3.partition()
        .size([2 * Math.PI, radius * radius]);

    // Nous dessinons un cercle complet pour que la détection de sortie du parent g (chart-container) ne se fasse pas au milieu du graphique
    vis.append("circle")
        .attr("r", radius)
        .style("opacity", 0);

    var root = d3.hierarchy(json)
        .sum(function (d) {
            return d.amount;
        })
        .sort(function (a, b) {
            if (a.depth === 1) { // Premier niveau trié par montant
                return b.value - a.value;
            } else { // Sinon trié par nom pour conserver la cohérence des patterns (depth = 2) et changement de saturation (depth = 3)
                return b.data.name.localeCompare(a.data.name) * -1;
            }
        });

    // Par efficacité, suppression des valeurs trop petites
    var nodes = partition(root).descendants()
        .filter(function (d) {
            return (d.x1 - d.x0 > 0.005); // 0.005 radians = 0.29 degrees
        });

    var path = vis.selectAll("path")
        .data(nodes)
        .enter().append("path")
        .attr("display", function (d) {
            return d.depth ? null : "none";
        })
        .attr("d", arc)
        .style("fill", function (d) {
            return getFillValue(d);
        })
        .on("mouseover", mouseover);

    // Ajout d'un mouseleave sur le parent g
    d3.select("#chart-container").on("mouseleave", mouseleave);
}

function getFillValue(d) {
    // Couleur principale de l'entreprise
    if (d.depth === 1) {
        return positionColors[d.data.position];
    }
    // Déclinaison de la couleur principale en patterns
    if (d.depth === 2) {
        let parentColor = positionColors[d.parent.data.position];
        // Pour obtenir un positionnement des patterns identiques sur le cercle, on fait varier la valeur rotate en fonction de la position de d
        // On considère que 1 radian (unité des arcs dans D3) = 57.29 degrée
        // Notre calcul permet d'avoir toujours le même pattern relativement au centre du cercle au milieu de l'arc dessiné
        let rotateValue = (d.x0 + d.x1) / 2 * 57.29;
        let patternId = d.data.name + d.parent.data.position;

        if (d.data.name == "Avantage") {
            let pattern = d3.select("#defs")
                .append("pattern")
                .attr("id", patternId)
                .attr("width", "8")
                .attr("height", "8")
                .attr("patternUnits", "userSpaceOnUse")
                .attr("patternTransform", "rotate(" + rotateValue + ")");
            pattern.append("rect")
                .attr("width", "4")
                .attr("height", "8")
                .attr("fill", parentColor);
            return "url(#" + patternId + ")";
        } else if (d.data.name == "Convention") {
            let pattern = d3.select("#defs")
                .append("pattern")
                .attr("id", patternId)
                .attr("width", "10")
                .attr("height", "10")
                .attr("patternUnits", "userSpaceOnUse");
            pattern.append("circle")
                .attr("r", "5")
                .attr("fill", parentColor);
            return "url(#" + patternId + ")";
        } else if (d.data.name == "Rémunération") {
            let pattern = d3.select("#defs")
                .append("pattern")
                .attr("id", patternId)
                .attr("width", "8")
                .attr("height", "8")
                .attr("patternUnits", "userSpaceOnUse")
                .attr("patternTransform", "rotate(" + (rotateValue + 90) + ")");
            pattern.append("rect")
                .attr("width", "1")
                .attr("height", "8")
                .attr("fill", parentColor);
            return "url(#" + patternId + ")";
        }
    }
    // Déclinaison de la couleur principale par rotation de Hue (saturation)
    if (d.depth === 3) {
        let parentColor = d3.hsl(positionColors[d.parent.parent.data.position]);
        parentColor.h += rotationColors[d.data.name];
        return parentColor + "";
    }

    return "";
}

function mouseover(d) {
    d3.select("#amount")
        .text(d.value.toLocaleString('fr-FR', { minimumFractionDigits: 0, style: 'currency', currency: 'EUR' }));

    var sequenceArray = d.ancestors().reverse();
    sequenceArray.shift(); // suppression de la racine

    d3.select("#category-amount")
        .text("");
    d3.select("#type-amount")
        .text("");

    sequenceArray.forEach(d => {
        if (d.depth === 1) {
            d3.select("#entreprise")
                .text(d.data.name);
        } else if (d.depth === 2) {
            d3.select("#type-amount")
                .text(d.data.name);
        } else if (d.depth === 3) {
            let text = d.data.name
                .replace("Académies, Fondation, sociétés savantes, organismes de conseils", "Académies, Fondation, ...")
                .replace("Personnes morales assurant la formation initiale ou continue des professionnels de santé", "Personnes morales assurant ...");
            d3.select("#category-amount")
                .text(text);
        }
    });

    // On grise tous les segments
    d3.selectAll("path")
        .style("opacity", 0.3);

    // Ensuite on met en valeur uniquement ceux qui sont ancêtres de la sélection
    vis.selectAll("path")
        .filter(function (node) {
            return (sequenceArray.indexOf(node) >= 0);
        })
        .style("opacity", 1);
}

function mouseleave(d) {
    // On désactive la fonction mouseover le temps de la transition
    d3.selectAll("path").on("mouseover", null);

    // Transition pour revenir à l'état d'origine et on remet le mouseover
    d3.selectAll("path")
        .transition()
        .duration(1000)
        .style("opacity", 1)
        .on("end", function () {
            d3.select(this).on("mouseover", mouseover);
        });
}

// Version minimaliste du Sunburst Chart
function buildMinChart(json) {
    const widthMin = 500, heightMin = 500, radiusMin = 250;

    const visMin = d3.select("#min-chart").append("svg")
        .attr("width", widthMin)
        .attr("height", heightMin)
        .append("g")
        .attr("transform", "translate(" + widthMin / 2 + "," + heightMin / 2 + ")");

    var arcMin = d3.arc()
        .startAngle(function (d) {
            return d.x0;
        })
        .endAngle(function (d) {
            return d.x1;
        })
        .innerRadius(function (d) {
            return Math.sqrt(d.y0);
        })
        .outerRadius(function (d) {
            return Math.sqrt(d.y1);
        });

    var root = d3.hierarchy(json)
        .sum(function (d) {
            return d.amount;
        });

    var partitionMin = d3.partition()
        .size([2 * Math.PI, radiusMin * radiusMin]);

    visMin.selectAll("path")
        .data(partitionMin(root).descendants())
        .enter().append("path")
        .attr("display", function (d) {
            return d.depth ? null : "none";
        })
        .attr("d", arcMin);
}

// From flat data to Company -> Type -> Category
function buildHierarchy(data) {
    var bigOnes = [
        { "displayName": "Sanofi", "matchingNames": ["sanofi"] },
        { "displayName": "Astrazeneca", "matchingNames": ["astrazeneca"] },
        { "displayName": "Novartis", "matchingNames": ["novartis"] },
        { "displayName": "Bristol-Myers Squibb", "matchingNames": ["bristol"] },
        { "displayName": "AbbVie", "matchingNames": ["abbvie"] },
        { "displayName": "Roche", "matchingNames": ["roche"] },
        { "displayName": "Merck Sharp and Dohme", "matchingNames": ["MSD France", "merck"] },
        { "displayName": "Microport CRM", "matchingNames": ["Microport"] },
        { "displayName": "GlaxoSmithKline", "matchingNames": ["glaxosmithkline"] }
        //{"displayName": "Lilly France", "matchingNames": ["lilly"]},
        //{"displayName": "Celgene", "matchingNames": ["celgene"]},
        //{"displayName": "Bayer", "matchingNames": ["bayer"]},
        //{"displayName": "Pfizer", "matchingNames": ["pfizer"]},
        //{"displayName": "Janssen-Cilag", "matchingNames": ["janssen"], "headName": "Johnson & Johnson"},
        //{"displayName": "Les laboratoires Servier", "matchingNames": ["servier"]},
        //{"displayName": "Applied Molecular Genetics", "matchingNames": ["amgen"]},
        //{"displayName": "Ipsen", "matchingNames": ["ipsen"]},
        //{"displayName": "Aurobindo Pharma", "matchingNames": ["arrow"]}
        //{"displayName": "Boehringer Ingelheim France", "matchingNames": ["boehringer"]},
        //{"displayName": "Gilead Sciences", "matchingNames": ["gilead"]},
        //{"displayName": "Groupe Pierre Fabre", "matchingNames": ["fabre"]},
        //{"displayName": "Medtronic France", "matchingNames": ["medtronic"]},
        //{"displayName": "Autres"}
    ];

    var root = { "name": "root", "children": [] };

    // Ajout des entreprises dont les montants sont les plus importants
    for (let i = 0; i < bigOnes.length; ++i) {
        root.children.push({
            "name": bigOnes[i].displayName, "matchingNames": bigOnes[i].matchingNames, "position": i, "children": [
                { "name": "Avantage", "children": [] }, {
                    "name": "Convention",
                    "children": []
                }, { "name": "Rémunération", "children": [] }
            ]
        });
    }

    for (let i = 0; i < data.length; ++i) {
        data[i].avant_montant = +data[i].avant_montant;
        data[i].conv_montant = +data[i].conv_montant;
        data[i].remu_montant = +data[i].remu_montant;

        let parentNode = getParentNode(root, data[i]);
        if (parentNode === undefined) {
            continue;
        }

        // Avantage
        if (data[i].avant_montant !== "" && data[i].avant_montant !== 0) {
            let foundCategory = undefined;
            for (var iCat = 0; iCat < parentNode.children[0].children.length; ++iCat) {
                if (parentNode.children[0].children[iCat].name === data[i].categorie) {
                    foundCategory = parentNode.children[0].children[iCat];
                }
            }
            if (foundCategory === undefined) {
                parentNode.children[0].children.push({ "name": data[i].categorie, "amount": 0 });
                foundCategory = parentNode.children[0].children[parentNode.children[0].children.length - 1];
            }
            foundCategory.amount = foundCategory.amount + data[i].avant_montant;
        }

        // Convention
        if (data[i].conv_montant !== "" && data[i].conv_montant !== 0) {
            let foundCategory = undefined;
            for (var iCat = 0; iCat < parentNode.children[1].children.length; ++iCat) {
                if (parentNode.children[1].children[iCat].name === data[i].categorie) {
                    foundCategory = parentNode.children[1].children[iCat];
                }
            }
            if (foundCategory === undefined) {
                parentNode.children[1].children.push({ "name": data[i].categorie, "amount": 0 });
                foundCategory = parentNode.children[1].children[parentNode.children[1].children.length - 1];
            }
            foundCategory.amount = foundCategory.amount + data[i].conv_montant;
        }

        // Rémunération
        if (data[i].remu_montant !== "" && data[i].remu_montant !== 0) {
            let foundCategory = undefined;
            for (var iCat = 0; iCat < parentNode.children[2].children.length; ++iCat) {
                if (parentNode.children[2].children[iCat].name === data[i].categorie) {
                    foundCategory = parentNode.children[2].children[iCat];
                }
            }
            if (foundCategory === undefined) {
                parentNode.children[2].children.push({ "name": data[i].categorie, "amount": 0 });
                foundCategory = parentNode.children[2].children[parentNode.children[2].children.length - 1];
            }
            foundCategory.amount = foundCategory.amount + data[i].remu_montant;
        }
    }

    return root;
}

function getParentNode(root, current) {
    let parentNode = undefined;
    for (let j = 0; j < root.children.length; ++j) {
        for (let k = 0; k < root.children[j].matchingNames.length; ++k) {
            if (current.denomination.search(new RegExp(root.children[j].matchingNames[k], "i")) >= 0) {
                return root.children[j];
            }
        }
    }
    return undefined;
    //return root.children[root.children.length - 1]; // Autres
}

function getSumCurrent(current) {
    return 0 + (current.avant_montant !== "" ? current.avant_montant : 0) + (current.conv_montant !== "" ? current.conv_montant : 0) + (current.remu_montant !== "" ? current.remu_montant : 0);
}
