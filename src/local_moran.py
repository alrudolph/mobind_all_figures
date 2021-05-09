import pandas as pd
import libpysal as lps
from esda import Moran_Local
from tqdm import tqdm

def false_discovery_rate(arr, sig):
    df = pd.DataFrame(arr, columns=["p"]).sort_values("p")
    df["i"] = np.arange(1, len(arr) + 1) * sig / len(arr)
    df["sig"] = (df["p"] < df["i"] + 0)
    return list(df.sort_index()["sig"])

def bonferroni(arr, sig):
    return list(np.array(arr) < sig / len(arr))

def filter_quadrants(arr):
    return [(a if a < 3 else 0) for a in arr]

def combine(sim, fdr, bon):
    return [
        (b + 4 if b != 0 else (f + 2 if f != 0 else s))
        for b, f, s in zip(bon, fdr, sim)
    ]

def n_permutations(df):
    return 999 # default value

def moran_quadrants(col, W, alpha, which):
    local_moran = Moran_Local(col, W, geoda_quads=True, permutations=n_permutations(col))

    ps = local_moran.p_sim
    qs = filter_quadrants(local_moran.q)

    if which == "fdr":
        f = false_discovery_rate(ps, alpha)
    elif which == "sim":
        f = [p < alpha for p in ps]
    elif which == "bon":
        f = bonferroni(ps, alpha)
    elif which == "all":
        fdr = false_discovery_rate(ps, alpha)
        bon = bonferroni(ps, alpha)
        sim = [p < alpha for p in ps]

        qs = combine(
            qs * np.array(sim),
            qs * np.array(fdr),
            qs * np.array(bon)
        )
        f = sim
    else:
        raise 'Valid p-value evaluations: "bon", "fdr", or "sim"'
        
    return list(qs * np.array(f))

def local_moran(grouped, cols, date_col: str, group_col: str, map_data: pd.DataFrame, map_group_col: str, limit=None, sig=0.05, which="fdr"):
    if not isinstance(cols, list):
        cols = [cols]

    output = pd.DataFrame(
        columns=[*cols, 'date']
    )

    W = lps.weights.Queen(map_data["geometry"])
    W.transform = 'r'

    for i, (date, group) in tqdm(enumerate(grouped), total=min(limit, len(grouped)) if limit else len(grouped)):

        if limit and i == limit:
            break 

        if group.isnull().values.any():
            continue

        # it's important we merge on map so grouping order is consistent
        ordered = pd.merge(
            map_data, group, left_on=map_group_col, right_on=group_col, how='left'
        )

        row = {
            'date': [min(group['date_start']), max(group['date_end'])] if ('date_start' in ordered.columns) else min(group[date_col])
        }

        for j, col in enumerate(cols):
            mq = moran_quadrants(ordered[col], W, sig, which=which)
            row[col] = mq
            
        output = output.append(row, ignore_index=True)

    return (output, map_data[map_group_col])


import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import geopandas as gpd
import math

class StackedChoropleth:
    
    def __init__(self, gpd_maps, export, cols, main_folder, region, desc="", titles='', plot_dim=None, fips_order=None):

        if not isinstance(gpd_maps, list):
            gpd_maps = [gpd_maps]

        if not isinstance(cols, list):
            cols = [cols]

        if not isinstance(titles, list):
            titles = [titles]

        if not isinstance(region, list):
            region = [region]

        if len(gpd_maps) > 1 and len(cols) > 1:
            raise "Can only specify multiple maps OR columns."

        plot_dim = (2, 1)    

        if plot_dim[0] * plot_dim[1] != len(cols) * len(gpd_maps) and plot_dim[0] * plot_dim[1] != len(cols) * len(gpd_maps) + 1:
            raise "Dimension does not match specified colums."

        self._gpd_maps = gpd_maps
        self._cols = cols
        self._export = export
        self._plot_dim = plot_dim
        self._desc = desc if desc else "-".join(cols)
        self._titles = titles if titles else cols
        self._region = region

        self._fips_order = fips_order
        self._state_col = "STATEFP"

        main_region_folder = f'{main_folder}/{"-".join(region)}'
        self.count_subfolder = f'{main_region_folder}/count/'
        self.recent_subfolder = f'{main_region_folder}/recent/'
        self.both_subfolder = f'{main_region_folder}/combined/'
        self.scatter_subfolder = f'{main_region_folder}/scatter/'

        self._count_labels = ["Low Count", "High Count"]
        self._recent_labels = ["1/1/2020", "12/31/2020"]
        self._count_title = "Number of Weeks"
        self._recent_title = "Recency"

        self._plot_title_size = 14
        self._font_size = 12

        self._figsize = (10, 8)

    #
    #   MAPPING OPTIONS
    #

    def plot_count(self):
        if '_collapse_count_hot' not in locals():
            self.__collapse_count()

        self.__create_choropleth(
            self._collapse_count_hot, 
            self._collapse_count_cold, 
            typ=self._count_title, 
            labels=self._count_labels, # Start date to End Date
            figsize=self._figsize,
            output_folder=self.count_subfolder
        )

    def plot_recent(self):
        if '_collapse_recent_hot' not in locals():
            self.__collapse_recent()
        
        self.__create_choropleth(
            self._collapse_recent_hot, 
            self._collapse_recent_cold, 
            typ=self._recent_title, 
            labels=self._recent_labels, # maybe actual values?
            figsize=self._figsize,
            output_folder=self.recent_subfolder
        )

    def plot_both(self):
        if '_collapse_count_combined' not in locals():
            self.__collapse_count_combined()

        if '_collapse_recent_hot' not in locals():
            self.__collapse_recent()
        
        hots = self._collapse_recent_hot
        colds = self._collapse_recent_cold

        fig, axes = plt.subplots(
            self._plot_dim[0], 
            self._plot_dim[1], 
            figsize=self._figsize,
        )
        fig.tight_layout()
        #plt.subplots_adjust(hspace=0.1)
    

        # Since I'm restricting maps OR col to be a single item this is really a single loop:
        for map_idx, gpd_map in enumerate(self._gpd_maps):

            for j, col in enumerate(self._cols):

                map_copy = gpd_map.copy()

                map_copy["geometry"] = [
                    col.centroid.buffer(10000 + 11000/40 * count)
                    for col, count in zip(map_copy['geometry'], self._collapse_count_combined[map_idx + j])
                ]
                # Row wise:
                ax = axes[map_idx + j]

                if self._region and self._region[map_idx] == "usa":
                    ax.set_xlim([-0.235e7, 0.22e7])
                    ax.set_ylim([-1.75e6, 1.45e6])
                elif self._region[map_idx] == "ca":
                    ax.set_xlim([-2.60e6, -1.563e6])
                    ax.set_ylim([-1e6, 0.65e6])
                elif self._region[map_idx] == "fl":
                    ax.set_xlim([0.730e6, 1.570e6])
                    ax.set_ylim([-1.700e6, -0.950e6])
                elif self._region[map_idx] == "ny":
                    ax.set_xlim([1.1e6, 2.0e6])
                    ax.set_ylim([.03e6, 0.9e6])
                elif self._region[map_idx] == "tx":
                    ax.set_xlim([-0.990e6, 0.24e6])
                    ax.set_ylim([-1.75e6, -0.380e6])

                # if adding new state start with looking at state bounds:
                # print(gpd_map.total_bounds)

                norm = colors.Normalize(vmin=1, vmax=max([*hots[map_idx + j], *colds[map_idx + j]]))

                if len(self._titles) > 1:
                    ax.set_title(self._titles[map_idx + j], fontsize=self._plot_title_size)
                # utility.hide_axis(ax)

                #self.__show_country_outline(ax, gpd_map)
                self.__show_state_outline(ax, gpd_map)
                self.__create_choropleth_map(hots[map_idx + j], ax, map_copy, self.__get_pallete("Reds"), norm)
                self.__create_choropleth_map(colds[map_idx + j], ax, map_copy, self.__get_pallete("Blues"), norm)
                
        if len(self._titles) == 1:
            fig.suptitle(self._titles[0], fontsize=self._plot_title_size)

        self.__create_choropleth_legend_horiz(fig, self._recent_title, self._recent_labels)
        self.__create_choropleth_legend_circles(fig, self._count_labels)
        # utility.save_plot(self._desc)

    #
    #   MAP FEATURES
    #

    def __create_choropleth(self, hots, colds, typ, labels, figsize, output_folder):
        fig, axes = plt.subplots(
            self._plot_dim[0], 
            self._plot_dim[1], 
            figsize=figsize
        )
        fig.tight_layout()
        # Since I'm restricting maps OR col to be a single item this is really a single loop:
        for map_idx, gpd_map in enumerate(self._gpd_maps):

            for j, col in enumerate(self._cols):
                # Row wise:
                #ax = utility.get_axis(axes, (map_idx + j) // self._plot_dim[0], (map_idx + j) % self._plot_dim[1])

                ax = axes[i]

                norm = colors.Normalize(vmin=0.5, vmax=max([*hots[map_idx + j], *colds[map_idx + j]]))
                ax.set_title(self._titles[map_idx + j], fontsize=self._plot_title_size)
                ax.set_axis_off()

                #self.__show_country_outline(ax, gpd_map)
                self.__show_state_outline(ax, gpd_map)
                self.__create_choropleth_map(hots[map_idx + j], ax, gpd_map, self.__get_pallete("Reds"), norm)
                self.__create_choropleth_map(colds[map_idx + j], ax, gpd_map, self.__get_pallete("Blues"), norm)
        
        self.__create_choropleth_legend_horiz(fig, typ, labels)

        # utility.save_plot(self._desc)

    def __create_choropleth_map(self, data, ax, gpd_map, palette, norm, **kwargs):
        gpd_map \
            .assign(cl = data) \
            .plot(column='cl', k=2, cmap=palette, edgecolor='black', ax=ax, linewidth=0, norm=norm, **kwargs)

    def __create_choropleth_legend_horiz(self, fig, typ, labels):
        #
        #   THIS NEEDS TO BE LOOKED INTO
        #
        start_x = 0.12 
        start_y = -0.07 # 0 for non tight
        
        newa = fig.add_axes([start_x, start_y + 0.04, 0.3, 0.03])
        sm = plt.cm.ScalarMappable(cmap="Reds")
        cb = plt.colorbar(sm, cax=newa, ticks=[], orientation="horizontal")
        cb.ax.tick_params(size=0)
        cb.ax.set_xlabel(typ, fontsize=self._font_size)
        cb.ax.xaxis.set_label_position("top")
        cb.ax.set_ylabel("Hot Spots", rotation=0, labelpad=35, y=0.125, fontsize=self._font_size)

        newa = fig.add_axes([start_x, start_y + 0, 0.3, 0.03])
        sm = plt.cm.ScalarMappable(cmap="Blues")
        cb = plt.colorbar(sm, cax=newa, ticks=[0.1, 0.875], orientation="horizontal")
        cb.ax.set_xticklabels(labels, fontsize=self._font_size)
        cb.ax.tick_params(size=0)
        cb.ax.set_ylabel("Cold Spots", rotation=0, labelpad=35, y=0.125, fontsize=self._font_size)

    def __create_choropleth_legend_vert(self, fig, typ, labels):
        #
        #   THIS NEEDS TO BE ADDED
        #
        return 1

    def __create_choropleth_legend_circles(self, fig, labels):
        start_y = -0.06 # 0 for non tight
        start_x = .58 # 0.5 for non tight
        
        if self._plot_dim[0] == 1 and self._plot_dim[1] == 1:
            point_scale = lambda i: (15000 + 55000/40 * i) / 2500
        else:
            point_scale = lambda i: (10000 + 11000/40 * i) / 2500

        newa = fig.add_axes([start_x, start_y + 0.07, 0.3, 0.03])
        points = [1, 5, 10, 20, 52]
        point_lgd = [plt.scatter([],[], s=point_scale(i), marker='o', color='k') for i in points]
        newa.legend(point_lgd, points, frameon=False, title=self._count_title, ncol=5, handlelength=0.1)
        newa.set_axis_off()


    #
    #   UTILITY MAPPING FUNCTIONS
    #

    def __show_country_outline(self, ax, gpd_map):
        gpd_map \
            .assign(dissolve = 0) \
            .dissolve(by="dissolve") \
            .plot(color="#FFFFFF00", ax=ax, edgecolor='black', linewidth=3)
        
    def __show_state_outline(self, ax, gpd_map):
        gpd_map \
            .dissolve(by=self._state_col) \
            .plot(color="#FFFFFF00", ax=ax, edgecolor='black', linewidth=1)

    def __get_pallete(self, which):
        import copy

        palette = copy.copy(plt.get_cmap(which))
        palette.set_under("#FFFFFF00", 0)
        return palette

    #
    #   CALCULATIONS:
    #
    def __collapse_count_combined(self):
        collapsed = {}

        for map_idx, gpd_map in enumerate(self._gpd_maps):
            for i, col in enumerate(self._cols):
                stacked = np.array([0] * gpd_map.shape[0])

                for week_arr in self._export[col]:
                    arr = [(a != 0) for a in self.__filter_state(gpd_map, week_arr)]
                    stacked = stacked + arr

                collapsed[map_idx + i] = stacked

        # set second title label to max number ??

        self._collapse_count_combined = collapsed

    def __collapse_count(self):
        collapsed_hot = {}
        collapsed_cold = {}

        for map_idx, gpd_map in enumerate(self._gpd_maps):
            for i, col in enumerate(self._cols):

                length = gpd_map.shape[0]

                stacked_hot = np.array([0] * length)
                stacked_cold = np.array([0] * length)

                for week_arr in self._export[col]:
                    arr_hot, arr_cold = list(zip(*[(a % 2 == 1, a != 0 and a % 2 == 0) for a in self.__filter_state(gpd_map, week_arr)]))
                    stacked_hot = stacked_hot + arr_hot
                    stacked_cold = stacked_cold + arr_cold

                collapsed_hot[map_idx + i] = stacked_hot
                collapsed_cold[map_idx + i] = stacked_cold

        # set second title label to max number ??

        self._collapse_count_hot = collapsed_hot
        self._collapse_count_cold = collapsed_cold
        
    def __collapse_recent(self):
        collapsed_hot = {}
        collapsed_cold = {}


        for map_idx, gpd_map in enumerate(self._gpd_maps):
            for i, col in enumerate(self._cols):

                length = gpd_map.shape[0]

                stacked_hot = [0] * length
                stacked_cold = [0] * length

                for week_num, week_arr in enumerate(self._export[col]):
                    arr_hot, arr_cold = list(zip(*[
                        ((week_num + 1) * (a % 2 == 1), (week_num + 1) * (a % 2 == 0 and a != 0)) 
                        for a in self.__filter_state(gpd_map, week_arr)
                    ]))
                    stacked_hot = [max(arr_hot[i], stacked_hot[i]) for i in range(length)]
                    stacked_cold = [max(arr_cold[i], stacked_cold[i]) for i in range(length)]

                # If hot and cold, pick most recent value ---- should be by cound I thinkk ahh
                for j, (h, c) in enumerate(zip(stacked_hot, stacked_cold)):
                    if h > 0 and c > 0:
                        if h > 0:
                            stacked_cold[j] = 0
                        else:
                            stacked_hot[j] = 0

                collapsed_hot[map_idx + i] = stacked_hot
                collapsed_cold[map_idx + i] = stacked_cold

        self._collapse_recent_hot = collapsed_hot
        self._collapse_recent_cold = collapsed_cold

    def __filter_state(self, map_data, arr):
    
        if self._fips_order is None:
            return arr
        
        data = pd.DataFrame()
        data["fips"] = self._fips_order
        data["val"] = arr

        return map_data.merge(data, how="left", on="fips")["val"]