    Because the computing resources allocated by GEE are insufficient to support the proposed algorithm, I have to divide it into several parts to run them separately. They can be opened using Notepad++ software. I take Fuzhou city (2000-2019) as an example to demonstrate the steps for running these algorithms:
(1)	Calculating the annual mean BCI (code: ¡°annual_Mean_BCI_calculation.js¡±). In the code, the period from 2010 to 2019 is set. The result from 2000 to 2009 can also be obtained by modifying the period and rerunning this code.
(2)	Calculating the trend of BCI (code: ¡°BCI_trend_calculation.js¡±);
(3)	Determining the period of NDVI (code: ¡°Monthly_mean_NDVI_calculation.js¡±);
(4)	Detecting the urbanized areas (code: ¡°Urbanized_area_detection.js¡±).

    If you have any questions, please contact: zhangyang2907@163.com.
