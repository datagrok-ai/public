import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  //Duplicates a string in Julia
  export async function juliaDup(s: string): Promise<string> {
    return await grok.functions.call('Samples:JuliaDup', { s });
  }

  export async function scatterPlotJulia(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:ScatterPlotJulia', { t, xColumnName, yColumnName, colorColumnName });
  }

  //Welch's t-test
  export async function tTestJulia(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
    return await grok.functions.call('Samples:TTestJulia', { data, x, y });
  }

  export async function octaveScript2(): Promise<any> {
    return await grok.functions.call('Samples:OctaveScript2', {});
  }

  //Displays contour plot
  export async function scottishHillsContourPlot(df: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Samples:ScottishHillsContourPlot', { df });
  }

  //Displays ridge plot
  export async function scottishHillsRidgePlot(df: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Samples:ScottishHillsRidgePlot', { df });
  }

  //Displays simple scatter plot
  export async function scottishHillsDemo(data: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Samples:ScottishHillsDemo', { data });
  }

  //Body Mass Index
  export async function bmi(height: number, weight: number): Promise<number> {
    return await grok.functions.call('Samples:BMI', { height, weight });
  }

  //Cell imaging segmentation based on Watershed algorithm
  export async function cellImagingSegmentation(file: DG.FileInfo): Promise<number> {
    return await grok.functions.call('Samples:CellImagingSegmentation', { file });
  }

  export async function datawigTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:DatawigTest', {});
  }

  //Gets EXIF data from PNG, JPEG, WEBP
  export async function exif(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Samples:EXIF', { file });
  }

  export async function glomEnvTestGlobal(): Promise<string> {
    return await grok.functions.call('Samples:GlomEnvTestGlobal', {});
  }

  export async function glomEnvTestInplace(): Promise<string> {
    return await grok.functions.call('Samples:GlomEnvTestInplace', {});
  }

  export async function glomEnvTest(): Promise<string> {
    return await grok.functions.call('Samples:GlomEnvTest', {});
  }

  //Image classification based of Efficientnet model
  export async function imageClassification(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Samples:ImageClassification', { file });
  }

  //Duplicates a string in Python
  export async function pythonDup(s: string): Promise<string> {
    return await grok.functions.call('Samples:PythonDup', { s });
  }

  //displays a list of locally installed python packages
  export async function pythonLibs(): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:PythonLibs', {});
  }

  //CWT (Morlet wavelet) scalogram plot
  export async function scalogramPython(data: DG.DataFrame, signal: DG.Column, sampleRate: number, w0: number, removeDc: boolean): Promise<any> {
    return await grok.functions.call('Samples:ScalogramPython', { data, signal, sampleRate, w0, removeDc });
  }

  export async function scatterPlotPython(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:ScatterPlotPython', { t, xColumnName, yColumnName, colorColumnName });
  }

  //Welch's t-test
  export async function tTestPython(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
    return await grok.functions.call('Samples:TTestPython', { data, x, y });
  }

  export async function ulmo(): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:Ulmo', {});
  }

  //One-way ANOVA
  export async function anova(table: DG.DataFrame, categories: DG.Column, variable: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:ANOVA', { table, categories, variable });
  }

  //Autocorrelation function (ACF) plots
  export async function acf(data: DG.DataFrame, columns: string[]): Promise<any> {
    return await grok.functions.call('Samples:ACF', { data, columns });
  }

  //Body Surface Area
  export async function bsa(height: number, weight: number): Promise<number> {
    return await grok.functions.call('Samples:BSA', { height, weight });
  }

  //Contour plot
  export async function contourPlot(data: DG.DataFrame, columns: string[]): Promise<any> {
    return await grok.functions.call('Samples:ContourPlot', { data, columns });
  }

  //Charting tool to create standard financial charts given a time series
  export async function financialChart(rates: DG.DataFrame, date: DG.Column, open: DG.Column, high: DG.Column, low: DG.Column, close: DG.Column, volume: DG.Column, adjusted: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:FinancialChart', { rates, date, open, high, low, close, volume, adjusted });
  }

  //Fit dose-response curve with a 4-parameter model with no constrain
  export async function fittingDRC(table: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:FittingDRC', { table, x, y });
  }

  //IIR Butterworth filter (high or low pass)
  export async function iirFilter(data: DG.DataFrame, signal: DG.Column, order: number, frequency: number, sampleRate: number, type: string): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:IIRFilter', { data, signal, order, frequency, sampleRate, type });
  }

  //Imputes (numerical) missing values using the kNN algorithm
  export async function knnR(data: DG.DataFrame, columns: string[], neighbours: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:KNNR', { data, columns, neighbours });
  }

  //Kolmogorov–Smirnov test (K–S test or KS test)
  export async function ksTest(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
    return await grok.functions.call('Samples:KSTest', { data, x, y });
  }

  //Linear Discriminant Analysis
  export async function lda(table: DG.DataFrame, predict: DG.Column, features: string[], perc: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:LDA', { table, predict, features, perc });
  }

  //Linear regression
  export async function linearRegression(data: DG.DataFrame, x: DG.Column, y: DG.Column, interceptZero: boolean): Promise<number> {
    return await grok.functions.call('Samples:LinearRegression', { data, x, y, interceptZero });
  }

  //Body Mass Index
  export async function listPackages(): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:ListPackages', {});
  }

  //Mixed Effects Model (LMER)
  export async function lmer(table: DG.DataFrame, features: string[], random: DG.Column, predict: DG.Column, perc: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:LMER', { table, features, random, predict, perc });
  }

  //MANOVA with 3 Dependent Variables ((x, y) ~ z)
  export async function manova(table: DG.DataFrame, variable1: DG.Column, variable2: DG.Column, variable3: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:MANOVA', { table, variable1, variable2, variable3 });
  }

  //Principal Component Analysis
  export async function pcaR(T: DG.DataFrame, columns: string[], numComp: number, center: boolean, scale: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:PCAR', { T, columns, numComp, center, scale });
  }

  //Partial Least Squares (PLS)
  export async function pls(table: DG.DataFrame, predict: DG.Column, features: string[], components: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:PLS', { table, predict, features, components });
  }

  //Predictive model based on Support vector machine (SVM)
  export async function predictiveModelSVM(datasetPredict: DG.DataFrame, predict: DG.Column, dataset: DG.DataFrame, columns: string[], fillMissing: boolean, perc: number): Promise<any> {
    return await grok.functions.call('Samples:PredictiveModelSVM', { datasetPredict, predict, dataset, columns, fillMissing, perc });
  }

  //Duplicates a string in R
  export async function rdup(s: string): Promise<string> {
    return await grok.functions.call('Samples:RDup', { s });
  }

  //Simple Spell Checking
  export async function renvSpellingExample(): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:RenvSpellingExample', {});
  }

  //CWT (Morlet wavelet) scalogram plot
  export async function scalogramR(data: DG.DataFrame, signal: DG.Column, sampleRate: number, octaves: number, voices: number, removeDc: boolean): Promise<any> {
    return await grok.functions.call('Samples:ScalogramR', { data, signal, sampleRate, octaves, voices, removeDc });
  }

  export async function scatterPlotR(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:ScatterPlotR', { t, xColumnName, yColumnName, colorColumnName });
  }

  //Sentiment classification of emotions and polarity
  export async function sentimentClassification(data: DG.DataFrame, col: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Samples:SentimentClassification', { data, col });
  }

  //Spectrogram plot
  export async function spectrogram(data: DG.DataFrame, signal: DG.Column, sampleRate: number, windowLength: number, timeStep: number, removeDc: boolean): Promise<any> {
    return await grok.functions.call('Samples:Spectrogram', { data, signal, sampleRate, windowLength, timeStep, removeDc });
  }

  //Surface plot
  export async function surfacePlot(t: DG.DataFrame, X: DG.Column, Y: DG.Column, Z: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:SurfacePlot', { t, X, Y, Z });
  }

  //Welch's t-test
  export async function tTestR(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
    return await grok.functions.call('Samples:TTestR', { data, x, y });
  }

  export async function ternaryPlot(t: DG.DataFrame, topColumnName: DG.Column, leftColumnName: DG.Column, rightColumnName: DG.Column, pointSize: number): Promise<any> {
    return await grok.functions.call('Samples:TernaryPlot', { t, topColumnName, leftColumnName, rightColumnName, pointSize });
  }

  //Time series decomposition
  export async function timeSeriesDecomposition(data: DG.DataFrame, dates: DG.Column, observations: DG.Column): Promise<any> {
    return await grok.functions.call('Samples:TimeSeriesDecomposition', { data, dates, observations });
  }

  //Time series forecasting using ARIMA
  export async function arimaForecasting(data: DG.DataFrame, dates: DG.Column, observations: DG.Column, P: number, D: number, Q: number, obsForecast: number): Promise<any> {
    return await grok.functions.call('Samples:ARIMAForecasting', { data, dates, observations, P, D, Q, obsForecast });
  }
}

export namespace queries {
  export async function accessOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:AccessOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function accessProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:AccessProducts', {});
  }

  export async function clickHouseStates(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseStates', {});
  }

  export async function clickHouseCountries(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseCountries', {});
  }

  export async function clickHouseProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseProducts', {});
  }

  export async function clickHouseEmployees(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseEmployees', {});
  }

  export async function clickHouseCustomers(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseCustomers', {});
  }

  export async function clickHouseOrderDetailsByQuantityProductNameCountry(quantity: number, productName: string, country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseOrderDetailsByQuantityProductNameCountry', { quantity, productName, country });
  }

  export async function clickHouseCustomersInCountry(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:ClickHouseCustomersInCountry', { country });
  }

  export async function mariaDBOrders(employeeid: number, shipvia: string, freight: number, shipcountry: string, shipcity: string, freightless1000: boolean, requireddate: any, orderdate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MariaDBOrders', { employeeid, shipvia, freight, shipcountry, shipcity, freightless1000, requireddate, orderdate });
  }

  export async function mariaDBProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MariaDBProducts', {});
  }

  export async function mssqlall(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLAll', {});
  }

  export async function mssqlbyInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByInt', { orderid });
  }

  export async function mssqlbyStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByStringPatternInt', { shipVia });
  }

  export async function mssqlbyDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByDouble', { freight });
  }

  export async function mssqlbyStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByStringChoices', { shipCountry });
  }

  export async function mssqlbyStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByStringPatternString', { shipCity });
  }

  export async function mssqlbyDatetime(requiredDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByDatetime', { requiredDate });
  }

  export async function mssqlbyStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLByStringPatternDatetime', { orderDate });
  }

  export async function mssqlorders(employeeId: number, shipVia: string, freight: string, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: string, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function mssqlproducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MSSQLProducts', { ProductID });
  }

  export async function mySQLOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MySQLOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function mySQLProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:MySQLProducts', {});
  }

  export async function oracleOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:OracleOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function oracleProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:OracleProducts', {});
  }

  export async function postgresStates(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresStates', {});
  }

  export async function postgresCountries(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresCountries', {});
  }

  export async function postgresProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresProducts', {});
  }

  export async function postgresEmployees(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresEmployees', {});
  }

  export async function postgresCustomers(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresCustomers', {});
  }

  export async function postgresOrderDetailsByQuantityProductNameCountry(quantity: number, productName: string, country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresOrderDetailsByQuantityProductNameCountry', { quantity, productName, country });
  }

  export async function postgresCustomersInCountry(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresCustomersInCountry', { country });
  }

  export async function postgresProductLookup(lookup: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresProductLookup', { lookup });
  }

  export async function postgresProductDetails(productName: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:PostgresProductDetails', { productName });
  }

  export async function ordersByEmployee(shipCountry: string, shipCity: string, customerId: string, employee: string): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:OrdersByEmployee', { shipCountry, shipCity, customerId, employee });
  }

  export async function domainUsage(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:DomainUsage', {});
  }

  export async function sqliteOrders(): Promise<DG.DataFrame> {
    return await grok.data.query('Samples:SQLiteOrders', {});
  }
}

export namespace funcs {

}
