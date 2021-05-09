import math
import matplotlib.pyplot as plt
import numpy as np

class SamplingBiasScatter:

    def __init__(self, export, cols, main_folder, region, desc="", titles='', labels='', plot_dim=None):
    
        if not isinstance(cols, list):
            cols = [cols]

        if not isinstance(titles, list):
            titles = [titles]

        if not plot_dim:
            if len(cols)  > 1:
                plot_dim = (2, math.ceil(len(cols)  / 2))
            else: 
                plot_dim = (1, 1)

        # plot dim can be 1 greater
        if plot_dim[0] * plot_dim[1] != len(cols) and plot_dim[0] * plot_dim[1] != len(cols) + 1:
            raise "Dimension does not match specified colums."

        self._cols = cols
        self._export = export
        self._plot_dim = plot_dim
        self._desc = desc if desc else "-".join(cols)
        self._titles = titles if titles else cols

        self._main_folder = f'{main_folder}/{region}/bias'

        self._xaxis = "pop"
        self._xlabel = "Census Population"
        self._label_col = "NAME"
        self._labels = labels

        self._font_size = 7
        self._title_size = 10

    def plot(self, figsize=(3, 2)):

        fig, axes = plt.subplots(self._plot_dim[0], self._plot_dim[1], figsize=figsize, sharex=True)


        for i, col in enumerate(self._cols):
            ax = axes

            self.__single_plot(self._labels[i], ax, col, self._titles[i])

        # utility save      
        fig.subplots_adjust(wspace=0.25, hspace=0.35)

 
    def __single_plot(self, ylabel, ax, col, title):

        ax.set_xlabel(self._xlabel, fontsize=self._font_size)
        ax.set_ylabel(ylabel, fontsize=self._font_size)

        df = self._export[[self._xaxis, col, self._label_col]].dropna()

        # plot points and line
        ax.scatter(df[self._xaxis], df[col], alpha = 0.5, s=60)
        xcoord, ycoord = self.__calc_line(col)
        ax.plot(xcoord, ycoord, color='k', linestyle='-', linewidth=1)

        ax.tick_params(axis='both', labelsize=self._font_size)
        ax.yaxis.get_offset_text().set_size(self._font_size * 2/3)
        ax.xaxis.get_offset_text().set_size(self._font_size * 2/3)

        for i, point in df.sort_values(by=col, ascending=False).iloc[0:5].iterrows():

            x_pos = point[self._xaxis] + (df[self._xaxis].max() - df[self._xaxis].min()) / 35
            point_ha = "left"

            # put these county labels to left of point
            if point[self._xaxis] > df[self._xaxis].max() * 0.9 or point["NAME"] in ["Maricopa, AZ", "Tarrant, TX"]:
                x_pos = point[self._xaxis] - (df[self._xaxis].max() - df[self._xaxis].min()) / 35
                point_ha = "right"

            # county name text
            ax.text(
                x_pos, point[col], 
                point[self._label_col], 
                verticalalignment='center',
                ha=point_ha,
                fontsize=self._font_size
            )
            
            r = self._export[[self._xaxis, col]].corr().loc[self._xaxis, col]

            # Correlation coefficient text
            ax.text(
                self._export[self._xaxis].min(), 
                max(ycoord) * .98, 
                f"Pearson R: {r:.2f}", 
                fontsize=self._font_size,
                fontweight="light"
            )

        ax.set_title(title, fontsize = self._title_size)


    def __calc_line(self, col):
        slope = self._export[col].sum() / self._export[self._xaxis].sum()
        xcoord = np.array([self._export[self._xaxis].min(), self._export[self._xaxis].max()])
        ycoord = slope * xcoord
        return (xcoord, ycoord)