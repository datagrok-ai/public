# MultiPlot

MultiPlot is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.

MultiPlot was developed in purpose to provide to user ability to easily create dashboards based on complex data which have common influencing parameter (usually presented as axis X on charts). Common axis X combined with compact form of presenting charts allow user quickly estimate data meaning.

Distribution of the vertical space among the plots is based on following principles:

1. There are fixed sized plots, size of which specified in pixels, for instance '50px';
2. There are proportional sized plots, size of which speicfied in percents, for instance '25%'
3. There is at least one plot, size of which equals all remained space.
4. Each plot could have title. Title could have default size specified in options or individual size.

Types of charts:
1. Scatterplot;
2. Line chart;
3. Bar chart;
4. Status markers (shape and color of markers are used for display data status)