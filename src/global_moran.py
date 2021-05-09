import pandas as pd
import libpysal as lps
from tqdm import tqdm
from esda import Moran

def global_moran(grouped, cols: str, date_col: str, group_col: str, map_data: pd.DataFrame, map_group_col: str, limit=None, sig=0.05):

    output = pd.DataFrame(
        columns=['date', *[i + "_est" for i in cols], *[i + "_p" for i in cols]]
    )
    W = lps.weights.Queen(map_data["geometry"])
    W.transform = 'r'

    for i, (date, group) in tqdm(enumerate(grouped), total=min(limit, len(grouped)) if limit else len(grouped)):

        if limit and i == limit:
            break 

        ordered = pd.merge(
            map_data, group, left_on=map_group_col, right_on=group_col, how='left'
        )

        row = {
            'date': [min(group['date_start']), max(group['date_end'])] if ('date_start' in ordered.columns) else min(group[date_col])
        }

        for j, col in enumerate(cols):
            values = ordered[[map_group_col, col, "geometry"]]
            mi = Moran(values[col], W, transformation='r', two_tailed=False, permutations=n_permutations(values))
            row[col + "_est"] = mi.I
            row[col + "_p"] = mi.p_sim
        output = output.append(row, ignore_index=True)

    return (output, map_data[map_group_col])

def n_permutations(df):
    return 999 # default value


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import math
import datetime as dt
from matplotlib.lines import Line2D

class GlobalMoranScatter:

    def __init__(self, export, cols, main_folder, region, desc="", titles='', plot_dim=None):
    
        if not isinstance(cols, list):
            cols = [cols]

        if not isinstance(titles, list):
            titles = [titles]

        
        plot_dim = (1, 2)

        # Plot dim can be 1 greater
        if plot_dim[0] * plot_dim[1] != len(cols) and plot_dim[0] * plot_dim[1] != len(cols) + 1:
            raise "Dimension does not match specified colums."

        self._cols = cols
        self._export = export
        self._plot_dim = plot_dim
        self._desc = desc if desc else "-".join(cols)
        self._titles = titles if titles else cols

        self._main_folder = f'{main_folder}/{region}/global'
        self._sig_level = 0.05

    def plot(self, figsize=(16, 16)):

        fig, axes = plt.subplots(self._plot_dim[0], self._plot_dim[1], figsize=figsize, sharex=True)

        for i, col in enumerate(self._cols):
            ax = axes[i]
        
            col_data = self._export[[col + "_est", col + "_p", "date"]] \
                .replace([np.inf, -np.inf], np.nan) \
            #    .dropna() # used when there is only one plot with missing values

            if col_data.empty:
                print(f"No values for {col}, moving on...")
                continue

            estimates = col_data[col + "_est"]
            sig = [(p <= self._sig_level) + 0 for p in col_data[col + "_p"]] # 1 = sig, 0 = insig

            start_dates = [i[0] for i in col_data["date"]]
           # dates = self._get_date_labels(start_dates) do own formatting

            self._single_scatter(ax, start_dates, estimates, sig, self._titles[i])

        # if doesn't fill bottom row set last plot off
        if self._plot_dim[0] * self._plot_dim[1] != len(self._cols):
            
            ax.set_axis_off()
            # colors = ['blue', 'red']
            colors = ['blue']
            lines = [Line2D([0], [0], color=c, marker='o') for c in colors]
            # labels = ['Significant', 'Not Significant']
            labels = ['Significant']
            ax.legend(lines, labels, loc='center')
            ax.set_title("", fontsize=12)
        
        fig.subplots_adjust(wspace=0.25, hspace=0.35)

    def _single_scatter(self, ax, dates, estimates, sig, title):
     #   xpoints = range(0, len(estimates))
        xpoints = dates
       # xpoints = mdates.drange(dates[0], dates[-1], dt.timedelta(days=1))

        # vertical line at march 11
        ax.axvline(x=dt.datetime.strptime("03-11-2020", "%m-%d-%Y"), ymin=0, ymax=1, color="gray")

        if len(np.unique(sig)) > 1:
            scatter = ax.scatter(xpoints, estimates, c=sig, zorder=10, cmap=colors.ListedColormap(['r', 'b']))
            plt.legend(handles=scatter.legend_elements()[0], labels=["Not Significant", "Significant"])
        else: 
            color = 'b' if sig[0] == 1 else 'r'
            label = 'Significant' if sig[0] == 1 else 'Not Significant'
            scatter = ax.scatter(xpoints, estimates, c=sig, zorder=10, cmap=colors.ListedColormap([color]))
            #plt.legend(handles=scatter.legend_elements()[0], labels=[label])

        # connect points
        ax.plot(xpoints, estimates, color="lightgray")

        # change x-axis ticks to date labels
#        ax.set_xticks(xpoints)
#        ax.set_xticklabels(dates, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_xticks([dt.datetime.strptime("03-01-2020", "%m-%d-%Y"), dt.datetime.strptime("06-01-2020", "%m-%d-%Y"), dt.datetime.strptime("09-01-2020", "%m-%d-%Y"), dt.datetime.strptime("12-01-2020", "%m-%d-%Y")])
        myFmt = mdates.DateFormatter('%b-%Y')
        ax.xaxis.set_major_formatter(myFmt)

        ax.set_title(title)
        #ax.set_xlabel("Week")
        ax.set_ylabel("Global Moran's I Estimate")

    def _get_date_labels(self, start_dates):
        months_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Noc", "Dec"]

        months = np.array([months_labels[i.month - 1] for i in start_dates])
        years = np.array([str(i.year) for i in start_dates])

        months[months == ['', *months[:-1]]] = ''
        years[years == ['', *years[:-1]]] = ''

        # fix first date overlap:
        if months[1] != '':
            months[0] = ''
            years[1] = years[0] if years[1] == '' else years[1]

        if years[1] != '':
            years[0] = ''

        # fix last date overlap:
        # TODO

        if len(np.unique(years)) > 2: # two years and ''s
            dates = [d + " " + y for d, y, in zip(months, years)]
        else:
            dates = months

        return dates