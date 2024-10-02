const fs = require('fs');
const path = require('path');

const exclude = ['meta', '.*tests']

const useLatest = (process.argv[2] === 'latest' || false);
const createFile = (process.argv[2] === 'create' || false);

const getDirectories = async source => (fs.readdirSync(source, {withFileTypes: true}))
    .filter(dirent => dirent.isDirectory())
    .map(dirent => dirent.name)

async function fileContainsContent(filename, str) {
    let exists = false;
    try {
        if (fs.existsSync(filename)) {
            const contents = await fs.promises.readFile(filename, 'utf-8');
            const matcher = new RegExp("(^|\\s)(" + str + ")($|\\s)", "mgi");
            exists = matcher.test(contents.toLowerCase());
        }
    } catch (err) {
        console.log(err);
    }
    return exists
}

async function getChangelog(filename, version) {
    const changelog = await fileContainsContent(filename, version);
    let summary = '';
    let changelogDate = '';
    if (changelog) {
        const content = fs.readFileSync(filename).toString();
        const summaryMatcher = new RegExp(`## ${version} \\((?:\\d{4}-\\d{2}-\\d{2}|WIP)\\)(?<text>[\\S\\s]*?)(?:$|### (?:Features|Bug|Fixed|Breaking)|##|\\*? ?Dependency)`, "g");
        let matches = Array.from(content.matchAll(summaryMatcher), x => x[1].trim());
        summary = matches[0];
        if (!summary) {
            const versionMatcher = new RegExp(`## ${version} \\((?:\\d{4}-\\d{2}-\\d{2}|WIP)\\)(?<text>[\\S\\s]*?)(?:$|## \\d+\\.\\d+\\.\\d+)`, "g");
            const versionMatch = Array.from(content.matchAll(versionMatcher), x => x[1].trim())[0];
            let linkRegex = /([\S ]*\[[\S ]+\]\(\S+\)[\S ]*)/g;
            matches = Array.from((versionMatch??"").matchAll(linkRegex), x => x[1].trim());
            summary = matches.join('\n');
            if (!summary) {
                const featuresRegex = /### Features:?(?<text>[\S\s]*?)(?:$|### (?:Features|Bug|Fixed|Breaking)|##|\*? ?Dependency)/g;
                matches = Array.from((versionMatch??"").matchAll(featuresRegex), x => x[1].trim());
                summary = matches.join('\n');
            }
            if (!summary) {
                const bugsRegex = /### Bug Fixes:?(?<text>[\S\s]*?)(?:$|### (?:Features|Bug|Fixed|Breaking)|##|\*? ?Dependency)/g;
                matches = Array.from((versionMatch??"").matchAll(bugsRegex), x => x[1].trim());
                summary = matches.join('\n');
            }
        }
        const dateMatcher = new RegExp(`## ${version} \\((?<date>\\d{4}-\\d{2}-\\d{2})\\)`, "g");
        matches = Array.from(content.matchAll(dateMatcher), x => x[1].trim());
        changelogDate = matches[0];
    }

    return [changelog, summary, changelogDate];
}

const rootDir = path.resolve('../');
const helpDirLoc = `${rootDir}/help/deploy/releases/plugins`;
const jsonTemplateLoc = `${helpDirLoc}/plugins.json`;

async function getPackages() {
    let packagesList = [];
    let jsonContent = [];
    if (useLatest) {
        try {
            jsonContent = JSON.parse(fs.readFileSync(jsonTemplateLoc, 'utf8'));
        } catch (err) {
            console.log(err);
        }
    }
    if (useLatest) {
        const packagesDirs = await getDirectories('../packages');
        for (const d of packagesDirs) {
            try {
                let packageJsonContent = JSON.parse(fs.readFileSync(`../packages/${d}/package.json`, 'utf8'));
                if (!jsonContent.find(p => (p.name === packageJsonContent.friendlyName || p.name === packageJsonContent.fullName) && p.version === packageJsonContent.version)) {
                    packagesList.push({
                        name: packageJsonContent.name,
                        version: packageJsonContent.version,
                        displayName: packageJsonContent.friendlyName || packageJsonContent.fullName,
                        dir: `packages/${d}`
                    })
                }
            } catch (e) {
                console.log(`WARN: ${d} is skipped because of occurred error: ${e}`)
            }
        }
    } else {
        const packagesSearch = await fetch('https://registry.npmjs.org/-/v1/search?text=@datagrok&size=1000').then(response => response.json()).then(result => result.objects.map(p => p.package.name));
        for (const p of packagesSearch) {
            packagesList.push({
                name: p
            })
        }
    }

    for (const p of packagesList.filter((p) => !exclude.find(value => {
        const matcher = new RegExp(value);
        return matcher.test(p.name);
    }))) {
        let latest = {};
        let response = {};
        let dir = p.dir;
        if (!p.displayName) {
            let status;
            latest = await fetch(`https://registry.npmjs.org/${p.name}/latest`).then(response => {
                status = response.ok;
                return response.json();
            });
            if (status && !latest.deprecated) {
                try {
                    dir = latest.repository.directory
                    p.displayName = latest.friendlyName || latest.fullName || p.name
                } catch (e) {
                    console.log(`WARN: ${p.name} is skipped because of occurred error: ${e}`)
                    continue;
                }
            } else {
                continue;
            }
        }
        if (p.version) {
            let status;
            const versionCheck = await fetch(`https://registry.npmjs.org/${p.name}/${p.version}`).then(response => {
                status = response.ok;
                return response.json();
            });
            if (status && !versionCheck.deprecated) {
                try {
                    dir = versionCheck.repository.directory;
                } catch (e) {
                    console.log(`WARN: ${p.name} is skipped because of occurred error: ${e}`);
                    continue;
                }
            } else {
                continue;
            }
        }
        response = await fetch(`https://registry.npmjs.org/${p.name}`).then(response => response.json());

        // const d = new Date();
        // d.setFullYear(d.getFullYear() - 1);


        if (p.version) {
            console.log(`Added plugin ${p.displayName} version ${p.version}`);
            const changelogInfo = await getChangelog(`${rootDir}/${dir}/CHANGELOG.md`, p.version);
            const changelog = changelogInfo[0];
            const summary = changelogInfo[1];
            const changelogDate = changelogInfo[2];
            jsonContent.push({
                name: p.displayName,
                date: response.time[p.version],
                version: p.version,
                summary: summary,
                dir: dir,
                changelog: changelog,
                changelogDate: changelogDate || response.time[p.version].split("T")[0],
                readme: fs.existsSync(`${rootDir}/${dir}/README.md`)
            });
        } else {
            console.log(`Added plugin ${p.displayName} with all published versions`)
            const regex = new RegExp('^[0-9]+\.[0-9]+\.[0-9]+$');
            for (const v of Object.entries(response.time)) {
                // if (regex.test(v[0]) && d > new Date(v[1])) {
                if (regex.test(v[0])) {
                    const changelogInfo = await getChangelog(`${rootDir}/${dir}/CHANGELOG.md`, v[0]);
                    const changelog = changelogInfo[0];
                    const summary = changelogInfo[1];
                    const changelogDate = changelogInfo[2];

                    jsonContent.push({
                        name: p.displayName,
                        date: v[1],
                        version: v[0],
                        summary: summary,
                        dir: dir,
                        changelog: changelog,
                        changelogDate: changelogDate || v[1].split("T")[0],
                        readme: fs.existsSync(`${rootDir}/${dir}/README.md`)
                    });
                }
            }
        }
    }

    jsonContent.sort(function (a, b) {
        const keyA = new Date(a.date), keyB = new Date(b.date);
        // Compare the 2 dates
        if (keyA > keyB) return -1;
        if (keyA < keyB) return 1;
        return 0;
    });

    return jsonContent;
}

if (createFile) {
    if (!fs.existsSync(jsonTemplateLoc)) {
        fs.writeFileSync(jsonTemplateLoc, JSON.stringify([], null, 2) + '\n');
    }
} else {
    getPackages().then((packages) => {
        fs.writeFileSync(jsonTemplateLoc, JSON.stringify(packages.slice(0, 1000), null, 2) + '\n');
        for (const p of packages.reduce(function (acc, cur) {
            return acc.filter(function (obj) {
                return obj.name !== cur.name
            }).concat([cur])
        }, [])) {
            const source = `${rootDir}/${p.dir}/CHANGELOG.md`
            const dest = `${helpDirLoc}/${p.name.replaceAll(' ', '_')}.md`
            if (fs.existsSync(source)) {
                fs.copyFile(source, dest, async function (err) {
                    if (err) {
                        console.log(err);
                    } else {
                        fs.readFile(dest, 'utf-8', async function (err, data) {
                            if (err) {
                                return console.log(err);
                            }
                            const exist = await fileContainsContent(dest, `title: ${p.name}`)
                            if (!exist) {
                                fs.open(dest, 'w+', function (err, fd) {
                                    if (err) throw err;
                                    const buffer = Buffer.from(`---\ntitle: ${p.name}\n---\n\n`);
                                    fs.write(fd, buffer, 0, buffer.length, 0, function (err) {
                                        if (err) throw err;
                                    });
                                    const dataBuffer = Buffer.from(data.replaceAll(/(\.\.\/)+help/g, '/help').replaceAll(/\.md/g, ''))
                                    fs.write(fd, dataBuffer, 0, dataBuffer.length, buffer.length, function (err) {
                                        if (err) throw err;
                                        fs.close(fd, function (err) {
                                            if (err) throw err;
                                        });
                                    });
                                });
                            }
                        });
                        console.log(`Added plugin ${p.name} changelog to help`)
                    }
                });
            }
        }
    })
}
