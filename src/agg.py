import libpysal as lps
from esda import Moran_Local, Moran
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
from tqdm import tqdm
import math

from typing import List

from datetime import datetime as dt

def get_year_week(date):
    """
    Gets the week number of date (starting on monday).
    Parameters
    ----------
    date : date
        Date to convert
    Returns
    -------
    tuple
        A tuple of (year, week_number)
    """
    week_num = date.isocalendar()[1]
    year = date.year

    # The week number can wrap to the next year, need to ensure year does as well.
    # Don't need to worry about wrapping other way
    if week_num == 1 and date.month == 12:
        year += 1

    if week_num == 53 and date.month == 1:
        year -= 1

    # Year comes first so dataframe is sorted chronologically
    return (year, week_num)

def _filters(n, coverage="usa", excl_non48=True):
    output = True

    if n == 11001: # DC
        output = False
    if excl_non48 and n >= 2000 and n <= 2999: # exclude AK
        output = False
    if excl_non48 and n >= 15001 and n <= 15009: # exclude HI
        output = False
    if n >= 60010: # territories and such
        output = False
    if n == 53055 or n == 25019: # ISLANDS WITH NO NEIGHBORS
        output = False
    if n == 51515: # Bedford County VA, code was changed in 2013
        output = False

    return output

def filter_data(data, coverage="usa", excl_non48=True):
    filters = [
        _filters(x["fips"], coverage, excl_non48) for _, x in data.iterrows()
    ]
    return data.loc[filters]

def filter_map(data, coverage="usa", excl_non48=True):
    # remove ? - I was doing something different for maps...

    filters = [
        _filters(x["fips"], coverage, excl_non48) for _, x in data.iterrows()
    ]
    return data.loc[filters]

def group_data(data, cols, date_col='', group_col='fips', by='week', date_format="%Y-%m-%d"):
    """
    Groups data by date to iterate over later. 

    Parameters
    ----------
    data : dataframe
        Must include a date and grouping column along with variables of interest
    cols : str | List[str]
        Column names for the variables of interest
    date_col : str
        Column name containing dates
    group_col : str
        Column name of grouping, i.e. fips code
    by: 'week' | 'day', optional
        How to group data, default: week.

    Returns
    -------
    dataframe
        Dataframe grouped by date
    """
    if not isinstance(cols, list):
        cols = [cols]

    if by != "county":
        # convert string dates to date objs
        dates = [dt.strptime(i, date_format).date() for i in data[date_col]]

    if by == 'week':
        year_week = [get_year_week(date) for date in dates]

        # aggregate the cols by mean, also get start (min) and end (max) dates
        agg_dict = dict(zip([*cols, "date_start", "date_end"],
                            [*np.repeat(np.nanmean, len(cols)), "min", "max"]))

        grouped = data \
            .assign(
                week_year=year_week,
                date_start=dates,
                date_end=dates
            ) \
            .groupby(["week_year", group_col]) \
            .agg(agg_dict) \
            .reset_index() \
            .groupby(["week_year"])

    elif by == 'day':
        grouped = data \
            .groupby([date_col])

        # In case there are multiple values for each day:
        # grouped = data \
        #    .assign(**{date_col: dates}) \
        #    .groupby([date_col, group_col])[cols].mean() \
        #    .groupby([date_col])

    elif by == 'county':
        grouped = data \
            .groupby([group_col]) \
            .agg(np.nanmedian) \
            .reset_index()

    return grouped


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib as mpl
import geopandas as gpd
import mapclassify
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

import platform
from pathlib import Path
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from matplotlib.lines import Line2D


def agvals(dframe, series, dt_col, id_col, time_level, loc_level, fips_filter=None):
    """
    Groups the data for further plotting by time and location. Filter by state or county fips. 
    
    Parameters:
    ---------- 
        dframe      (pd.DataFrame) - dataframe
        series      (str)          - which series to plot 
        dt_col      (str)          - datetime column 
        id_col      (str)          - id column 
        time_level  (str)          - either 'day' or 'week'
        loc_level   (str)          - either 'state' or 'county'
        fips_filter (str)          - state / county fips 2/5 digit id
    
    Returns: 
        grouped (pd.DataFrame)
        
    """
    
    # convert dt to datetime 
    dframe[dt_col] = pd.to_datetime(dframe[dt_col])
    
    # reindex df for easier grouping
    dframe.index = dframe[dt_col]
    
    # groupby statements
    # week AND county
    if (time_level=='week')&(loc_level=='county'):
        group_arr = [pd.Grouper(freq='W'), dframe[id_col], dframe["Region"]] if "Region" in dframe.columns else  [pd.Grouper(freq='W'), dframe[id_col]]
        grouped = ( dframe \
            .groupby(group_arr)[series] \
            .mean() \
            .reset_index())

    # week AND state
    elif (time_level=='week')&(loc_level=='state'):
        grouped = dframe \
            .groupby([pd.Grouper(freq='W'), dframe[id_col].str.slice(start=0,stop=2)])[series] \
            .mean() \
            .reset_index()
        
    # day and county
    # do nothing --> no grouping 
    elif (time_level=='day')&(loc_level=='county'):
        grouped = dframe
        grouped = grouped[[dt_col, id_col, series]]
    
    # day AND state 
    elif (time_level=='day')&(loc_level=='state'):
        grouped = dframe \
            .groupby([dframe[id_col].str.slice(start=0,stop=2)])[series] \
            .mean() \
            .reset_index()
        
    
    # filter
    if len(str(fips_filter))==2:
        grouped = grouped.loc[grouped[id_col].str.slice(start=0, stop=2)==fips_filter,]
    
    if len(str(fips_filter))==5: 
        grouped = grouped.loc[grouped[id_col]==fips_filter,]
    
    # rename cols for easier calling 
    grouped.columns = ['ds', 'fips', "Region", 'mobind'] if "Region" in dframe else ['ds', 'fips', 'mobind']
    
    # add state fips code 
    grouped['num_code'] = grouped.fips.str.slice(start=0, stop=2).astype(str)
    
    # calculate week number 
    grouped['w'] = grouped.ds.dt.week

  #  print(dframe.head(), grouped.head())

    return grouped

def plot_moby(dframe, path, title=None, figw=16, figh=6, kind='triple'):
    """
    Plots the data.

    Parameters:
    ----------- 
    dframe : dataframe
        pandas dataframe with data to be plotted 

    path: string 
        path where to store the jpeg 

    title: string 
        title for the figure 

    figw : int
        figure width

    figh : int 
        figure height

    kind : str
        Type of figure. One of 
        'triple' - time series plot 
        'calendar' - calendar heatmap showing the number of null values  
        'choro' - choropleth map showing the number of null values  

    Returns: 
        matplotlib figure 

    """

    fig, ax = plt.subplots(figsize=(figw, figh))

    g = sns.lineplot(data=dframe, x="ds", y="mobind", ax=ax, legend=True, lw=2)
    ca = sns.lineplot(data=dframe.loc[dframe.num_code=='06'], x="ds", y="mobind", 
                  lw=2, ax=ax, legend=True)
    sb = sns.lineplot(data=dframe.loc[dframe.fips=='06083'], x="ds", y="mobind", 
                  lw=2, ax=ax, legend=True)

    colors = ['tab:blue', 'tab:orange', 'tab:green']
    lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
    labels = ['USA', 'California', 'Santa Barbara county']

    ax.set(xlabel='', ylabel='')
    
    if title is not None:
        ax.set_title(title)

    plt.legend(lines, labels)

    if path is not None: 
        fig.savefig(path, bbox_inches='tight', dpi=150)
        plt.close(fig)

    return fig

# source: https://stackoverflow.com/questions/32485907/matplotlib-and-numpy-create-a-calendar-heatmap

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


DAYS = ['Sun.', 'Mon.', 'Tues.', 'Wed.', 'Thurs.', 'Fri.', 'Sat.']
MONTHS = ['Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'June', 'July', 'Aug.', 'Sept.', 'Oct.', 'Nov.', 'Dec.']


def date_heatmap(series, start=None, end=None, mean=False, ax=None, **kwargs):
    '''Plot a calendar heatmap given a datetime series.

    Arguments:
        series (pd.Series):
            A series of numeric values with a datetime index. Values occurring
            on the same day are combined by sum.
        start (Any):
            The first day to be considered in the plot. The value can be
            anything accepted by :func:`pandas.to_datetime`. The default is the
            earliest date in the data.
        end (Any):
            The last day to be considered in the plot. The value can be
            anything accepted by :func:`pandas.to_datetime`. The default is the
            latest date in the data.
        mean (bool):
            Combine values occurring on the same day by mean instead of sum.
        ax (matplotlib.Axes or None):
            The axes on which to draw the heatmap. The default is the current
            axes in the :module:`~matplotlib.pyplot` API.
        **kwargs:
            Forwarded to :meth:`~matplotlib.Axes.pcolormesh` for drawing the
            heatmap.

    Returns:
        matplotlib.collections.Axes:
            The axes on which the heatmap was drawn. This is set as the current
            axes in the `~matplotlib.pyplot` API.
    '''
    # Combine values occurring on the same day.
    dates = series.index.floor('D')
    group = series.groupby(dates)
    series = group.mean() if mean else group.sum()

    # Parse start/end, defaulting to the min/max of the index.
    start = pd.to_datetime(start or series.index.min())
    end = pd.to_datetime(end or series.index.max())

    # We use [start, end) as a half-open interval below.
    end += np.timedelta64(1, 'D')

    # Get the previous/following Sunday to start/end.
    # Pandas and numpy day-of-week conventions are Monday=0 and Sunday=6.
    start_sun = start - np.timedelta64((start.dayofweek + 1) % 7, 'D')
    end_sun = end + np.timedelta64(7 - end.dayofweek - 1, 'D')

    # Create the heatmap and track ticks.
    num_weeks = (end_sun - start_sun).days // 7
    heatmap = np.zeros((7, num_weeks))
    ticks = {}  # week number -> month name
    for week in range(num_weeks):
        for day in range(7):
            date = start_sun + np.timedelta64(7 * week + day, 'D')
            if date.day == 1:
                ticks[week] = MONTHS[date.month - 1]
            if date.dayofyear == 1:
                ticks[week] += f'\n{date.year}'
            if start <= date < end:
                heatmap[day, week] = series.get(date, 0)

    # Get the coordinates, offset by 0.5 to align the ticks.
    y = np.arange(8) - 0.5
    x = np.arange(num_weeks + 1) - 0.5

    # Plot the heatmap. Prefer pcolormesh over imshow so that the figure can be
    # vectorized when saved to a compatible format. We must invert the axis for
    # pcolormesh, but not for imshow, so that it reads top-bottom, left-right.
    ax = ax or plt.gca()
    mesh = ax.pcolormesh(x, y, heatmap, **kwargs)
    ax.invert_yaxis()

    # Set the ticks.
    ax.set_xticks(list(ticks.keys()))
    ax.set_xticklabels(list(ticks.values()))
    ax.set_yticks(np.arange(7))
    ax.set_yticklabels(DAYS)

    # Set the current image and axes in the pyplot API.
    plt.sca(ax)
    plt.sci(mesh)

    return ax


def date_heatmap_demo():
    '''An example for `date_heatmap`.

    Most of the sizes here are chosen arbitrarily to look nice with 1yr of
    data. You may need to fiddle with the numbers to look right on other data.
    '''
    # Get some data, a series of values with datetime index.
    data = np.random.randint(5, size=365)
    data = pd.Series(data)
    data.index = pd.date_range(start='2017-01-01', end='2017-12-31', freq='1D')

    # Create the figure. For the aspect ratio, one year is 7 days by 53 weeks.
    # We widen it further to account for the tick labels and color bar.
    figsize = plt.figaspect(7 / 56)
    fig = plt.figure(figsize=figsize)

    # Plot the heatmap with a color bar.
    ax = date_heatmap(data, edgecolor='black')
    plt.colorbar(ticks=range(5), pad=0.02)

    # Use a discrete color map with 5 colors (the data ranges from 0 to 4).
    # Extending the color limits by 0.5 aligns the ticks in the color bar.
    cmap = mpl.cm.get_cmap('Blues', 5)
    plt.set_cmap(cmap)
    plt.clim(-0.5, 4.5)

    # Force the cells to be square. If this is set, the size of the color bar
    # may look weird compared to the size of the heatmap. That can be corrected
    # by the aspect ratio of the figure or scale of the color bar.
    ax.set_aspect('equal')

    # Save to a file. For embedding in a LaTeX doc, consider the PDF backend.
    # http://sbillaudelle.de/2015/02/23/seamlessly-embedding-matplotlib-output-into-latex.html
    fig.savefig('heatmap.pdf', bbox_inches='tight')

    # The firgure must be explicitly closed if it was not shown.
    plt.close(fig)