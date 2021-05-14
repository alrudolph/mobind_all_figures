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
