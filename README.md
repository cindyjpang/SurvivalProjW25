# Survival Project  

## File Organization. 

### python > biostatProjData.ipynb. 
retrieves the original data files and performs merging for phenotype data. Exports to a .csv file for final data processing in R.  

### r > data. 
- **survDataCleaned.csv** contains the data relevant to running survival models.  
- **survProjData.csv** is the unprocessed phenotype and sample data. 

### r > src. 
- **models.R** has all the fitted full models and models from Stepwise AIC and BIC selection. To pull up all the models make sure you source from this file using `source(../src/models.R)`. 
- **CoxSnellResidualPlot.R** contains functions for plotting Cox-Snell residuals. 

### r > run.  
- **generateCSPlots.R** generates the Cox-Snell Residual Plots using the functions from `CoxSnellResidualPlot.R` and models from `models.R`. Plots are exported to r > figs. 
- **runModelSelection.Rmd** generates a pdf that shows the final models that were selected using stepwise AIC and BIC selection. See `runModelSelection.pdf` for results.  


