// String similarity as normalized distances in [0, 1] (0 = identical)

grok.shell.info(DG.StringUtils.levenshteinDistance('kitten', 'sitting')); // 0.4286
grok.shell.info(DG.StringUtils.jaroWinklerDistance('martha', 'marhta')); // 0.0389
