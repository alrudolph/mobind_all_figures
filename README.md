# Assessing COVID-induced changes in spatiotemporal structure of mobility in the United States in 2020: A multi-source analytical framework.

Evgeny Noi (noi@ucsb.edu), Alexander Rudolph, Somayeh Dodge

## Abstract

The COVID-19 pandemic resulted in profound changes in mobility patterns and altered travel behaviors locally and globally. As a result, movement metrics have widely been used by researchers and policy makers as indicators to study, model, and mitigate the impacts of the COVID-19 pandemic. However, the veracity and variability of these mobility metrics have not been studied. This paper provides a systematic review of mobility and social distancing metrics available to researchers during the pandemic in 2020 in the United States. Twenty-six indices across nine different sources are analyzed and assessed with respect to their spatial and temporal coverage as well as sample representativeness at the county-level. Finally global and local indicators of spatial association are computed to explore spatial and temporal heterogeneity in mobility patterns. The structure of underlying changes in mobility and social distancing is examined in different US counties and across different data sets. We argue that a single measure might not describe all aspects of mobility perfectly. A more comprehensive measure of mobility is required to model the complex effects of the epidemic. We urge that various limitations inherent in data sets from different providers should be acknowledged and these data must be used with caution.

```
project
|   counties.json
|   sg_data.csv
|   state_names.csv
|   state_regions.csv
|
|   figures1.ipynb
|   figures1.html
|
|___src
    |   agg.py
    |   global_moran.py
    |   local_moran.py
    |   sampling_bias.py
    |   violin.py
```

**figures1.ipynb** includes all of the code to generate figures used for the paper, also included in html format.

**figures1.html** inlcudes exported ipynb notebook for, for better viewing experience, than ipynb.

## Data

**counties.json** contains county geometries and data for our geopandas dataframe.

**sg_data.csv** contains the Safegraph daily county level data for this notebook.

**state_names.csv** is used to get the full state names from their abbreviations.

**state_regions** is used to classify states in the Northeast, West, South or Midwest.

## src/

### agg.py

Includes aggregation functions for the data

### global_moran.py

Includes functions to calculate global moran I test static and p-value, and plot those values

### local_moran.py

Includes functions to compute the local moran I test static, quadrants, etc... and make the bubble maps of the US

### sampling_bias.py

Includes function to make the plot of the sampling bias

### violin.py

Includes functions to scale the data and axis of the violin plots. The plots are in the notebook.
