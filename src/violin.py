
# modified from : http://maxberggren.se/2016/11/21/right-labels/

def legend_positions(dfs, ymin, ymax):
    """ Calculate position of labels to the right in plot... """
    positions = {}
    for i, df in enumerate(dfs):
        # grouped = dfs.groupby("ds", as_index=False).agg({ "mobind": np.nanmean }).dropna()["mobind"][-]
        vals = df.groupby("ds", as_index=False).agg({ "mobind": np.nanmean }).dropna()
        positions[i + 1] = np.array(vals["mobind"])[-1]

    def push():
        """
        ...by puting them to the last y value and
        pushing until no overlap
        """
        collisions = -1
        
        while(collisions != 0):
            collisions = 0

            for column1, value1 in positions.items():
                for column2, value2 in positions.items():

                    if column1 != column2:
                        dist = abs(value1-value2)
                        tot = abs(ymax - ymin)

                        # print(dist/tot)

                        if dist/tot < 0.1:
                            collisions += 1
                            if value1 < value2:
                                positions[column1] -= .01 * tot
                                positions[column2] += .01 * tot
                            else:
                                positions[column1] += .01 * tot
                                positions[column2] -= .01 * tot
                            return True
    while True:
        pushed = push()
        if not pushed:
            break

    print(positions)

    return positions


from sklearn.preprocessing import RobustScaler, StandardScaler
def scaleData(arr):
    a = StandardScaler(
        #with_mean=False, 
        # with_scaling=True,
        # unit_variance=True
    ).fit_transform(
        np.array(arr)[:, np.newaxis]
    )
    return np.concatenate(a, axis=0)

import numpy as np
from numpy import ma
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import FixedLocator, FuncFormatter

def set_outlier_cutoff(outlier_cutoff, minval, maxval):
    class MercatorLatitudeScale(mscale.ScaleBase):
        """
        Scales data in range -pi/2 to pi/2 (-90 to 90 degrees) using
        the system used to scale latitudes in a Mercator__ projection.

        The scale function:
        ln(tan(y) + sec(y))

        The inverse scale function:
        atan(sinh(y))

        Since the Mercator scale tends to infinity at +/- 90 degrees,
        there is user-defined threshold, above and below which nothing
        will be plotted.  This defaults to +/- 85 degrees.

        __ http://en.wikipedia.org/wiki/Mercator_projection
        """

        # The scale class must have a member ``name`` that defines the string used
        # to select the scale.  For example, ``gca().set_yscale("mercator")`` would
        # be used to select this scale.
        name = 'outlier_cutoff'

        def __init__(self, axis, *, thresh=np.deg2rad(85), **kwargs):
            """
            Any keyword arguments passed to ``set_xscale`` and ``set_yscale`` will
            be passed along to the scale's constructor.

            thresh: The degree above which to crop the data.
            """
            super().__init__(axis)
            self.thresh = thresh

        def get_transform(self):
            """
            Override this method to return a new instance that does the
            actual transformation of the data.

            The MercatorLatitudeTransform class is defined below as a
            nested class of this one.
            """
            return self.MercatorLatitudeTransform(self.thresh)

        def set_default_locators_and_formatters(self, axis):
            """
            Override to set up the locators and formatters to use with the
            scale.  This is only required if the scale requires custom
            locators and formatters.  Writing custom locators and
            formatters is rather outside the scope of this example, but
            there are many helpful examples in :mod:`.ticker`.

            In our case, the Mercator example uses a fixed locator from -90 to 90
            degrees and a custom formatter to convert the radians to degrees and
            put a degree symbol after the value.
            """
            #axis.set(major_locator=FixedLocator(np.array([-5, 0, 5, 10, 15, (maxval // 5 + 1) * 5])))

            # upper = max([
            #     (maxval // 5 + 1) * 5, 
            #     outlier_cutoff
            # ])

            # lower = min([
            #     math.floor(minval / 5 - 1) * 5, 
            #     -outlier_cutoff
            # ])

            # max/min val is next multiple of 10
            upper = [] if maxval < outlier_cutoff else [(maxval // 5 + 1) * 5 ]
            lower = [] if minval > -outlier_cutoff else [math.floor(minval / 5) * 5]

            r = np.array([*lower, -5, 0, 5, *upper])
            print("RANGE", r)

            axis.set(major_locator=FixedLocator(r))

            # if maxval > outlier_cutoff + 5:
            #     axis.set(major_locator=FixedLocator(np.array([
            #         *list(range(-5, outlier_cutoff + 5, 5)), (maxval // 5 + 1) * 5
            #     ])))
            # else:
            #     axis.set(major_locator=FixedLocator(
            #         np.array(range(
            #             (math.ceil(minval / 5) -1) * 5, 
            #             outlier_cutoff + 5 + (5 if maxval > outlier_cutoff else 0)
            #         , 5))
            #     ))

        def limit_range_for_scale(self, vmin, vmax, minpos):
            """
            Override to limit the bounds of the axis to the domain of the
            transform.  In the case of Mercator, the bounds should be
            limited to the threshold that was passed in.  Unlike the
            autoscaling provided by the tick locators, this range limiting
            will always be adhered to, whether the axis range is set
            manually, determined automatically or changed through panning
            and zooming.
            """
            return vmin, vmax

        class MercatorLatitudeTransform(mtransforms.Transform):
            # There are two value members that must be defined.
            # ``input_dims`` and ``output_dims`` specify number of input
            # dimensions and output dimensions to the transformation.
            # These are used by the transformation framework to do some
            # error checking and prevent incompatible transformations from
            # being connected together.  When defining transforms for a
            # scale, which are, by definition, separable and have only one
            # dimension, these members should always be set to 1.
            input_dims = output_dims = 1

            def __init__(self, thresh):
                mtransforms.Transform.__init__(self)
                self.thresh = thresh

            def transform_non_affine(self, a):
                """
                This transform takes a numpy array and returns a transformed copy.
                Since the range of the Mercator scale is limited by the
                user-specified threshold, the input array must be masked to
                contain only valid values.  Matplotlib will handle masked arrays
                and remove the out-of-range data from the plot.  However, the
                returned array *must* have the same shape as the input array, since
                these values need to remain synchronized with values in the other
                dimension.
                """

                upper_length = maxval - outlier_cutoff
                lower_length = -outlier_cutoff - minval
                total = maxval - minval

                def scale(o):
                    # dividing by max(upper_lengt, 5) so that outlier space doesn't get too big for just barely outliers
                    # multiplying by min(30, total) to keep the outlier space a minimum size

                    if o > outlier_cutoff:
                        return outlier_cutoff + (o - outlier_cutoff) / max(upper_length, 10) * min(10, total) / 5
                    elif o < -outlier_cutoff:
                        return -outlier_cutoff - (-outlier_cutoff - o) /  max(lower_length, 10) * min(10, total) / 5
                    return o

                return np.array([scale(o) for o in a ])

                # if maxval > outlier_cutoff + 5:
                #     lower_length = outlier_cutoff - minval + 5
                #     outlier_length = maxval - outlier_cutoff
                #     return np.array([
                #         o # ((o - minval) / lower_length * 4/5  
                #         if o < outlier_cutoff 
                #         else outlier_cutoff + (o - outlier_cutoff)/outlier_length * 5/24 * lower_length #/ 5#* lower_length / 5
                #         for o in a
                #     ])
                # else:
                #     return a

            def inverted(self):
                """
                Override this method so Matplotlib knows how to get the
                inverse transform for this transform.
                """
                return MercatorLatitudeScale.InvertedMercatorLatitudeTransform(
                    self.thresh)

        class InvertedMercatorLatitudeTransform(mtransforms.Transform):
            input_dims = output_dims = 1

            def __init__(self, thresh):
                mtransforms.Transform.__init__(self)
                self.thresh = thresh

            def transform_non_affine(self, a):
                return a

            def inverted(self):
                return MercatorLatitudeScale.MercatorLatitudeTransform(self.thresh)
    mscale.register_scale(MercatorLatitudeScale)