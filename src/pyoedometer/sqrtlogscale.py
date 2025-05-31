import warnings
import pdb
import sys
import numpy as np
from numpy import ma
import math
import matplotlib.scale as mscale
from matplotlib import rcParams
from matplotlib.transforms import Transform, IdentityTransform, nonsingular
from matplotlib.ticker import (NullFormatter, ScalarFormatter,
                               LogFormatterSciNotation, LogitFormatter)
from matplotlib.ticker import (NullLocator, LogLocator, AutoLocator,
                               SymmetricalLogLocator, LogitLocator, Locator)

try:
    from matplotlib.ticker import is_decade, decade_up, decade_down
except ImportError:
    from matplotlib.ticker import (
        _is_decade as is_decade,
        _decade_greater as decade_up,
        _decade_less   as decade_down
    )


import logging

import six
if six.PY3:
    long = int

class RootTransformBase(Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, nonpos):
        Transform.__init__(self)
        if nonpos == 'mask':
            self._fill_value = np.nan
        else:
            self._fill_value = 1e-300

    def transform_non_affine(self, a):
        with np.errstate(invalid="ignore"):
            a = np.where(a < 0, self._fill_value, a)
        return ma.power(a, 1./self.base)


class InvertedRootTransformBase(Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def transform_non_affine(self, a):
        return ma.power(a, self.base)


class SqrtTransform(RootTransformBase):
    base = 2

    def inverted(self):
        return InvertedSqrtTransform()


class InvertedSqrtTransform(InvertedRootTransformBase):
    base = 2

    def inverted(self):
        return SqrtTransform()


class CbrtTransform(RootTransformBase):
    base = 3

    def inverted(self):
        return InvertedCbrtTransform()


class InvertedCbrtTransform(InvertedRootTransformBase):
    base = 3

    def inverted(self):
        return CbrtTransform()


class SqrtLogTransform(Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    logbase = 10.
    rootbase = 2.

    def __init__(self, nonpos='mask', intersect=1., linscale=1.):
        Transform.__init__(self)
        if nonpos == 'mask':
            self._fill_value = np.nan
        else:
            self._fill_value = 1e-300
        self.intersect = intersect
        self.linscale = linscale

    def transform_non_affine(self, a):
        with np.errstate(invalid="ignore"):
            a = np.where(a < 0, self._fill_value, a)
            b = np.where(a <= 0, self._fill_value, a)

        return np.where(a <= self.intersect,
                        ma.power(a, 1./self.rootbase)/ma.power(self.intersect, 1./self.rootbase),
                        np.log(b/self.intersect, out=b)/np.log(self.logbase)+1)
                        #np.log10(b/self.intersect, out=b)+ma.power(self.intersect, 1./2.)/ma.power(self.intersect, 1./2.))

        # The division by sqrt(intersect) scales the sqrt section of the scale to
        # the length of approximately one decade on the log part of the scale

    def inverted(self):
        return InvertedSqrtLogTransform(intersect=self.intersect)


class InvertedSqrtLogTransform(Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    logbase = 10.
    rootbase = 2.

    def __init__(self, intersect=1.):
        Transform.__init__(self)
        self.intersect = intersect

    def transform_non_affine(self, a):
        return np.where(a < ma.power(self.intersect, 1./self.rootbase),
                        ma.power(a * np.power(self.intersect, 1./self.rootbase) , self.rootbase),
                        np.power(self.logbase, a - 1 + np.log(self.intersect)/np.log(self.logbase)))
                        #ma.power(a, 2.)*ma.power(self.intersect, 1./2.)/ma.power(self.intersect, 1./2.),
                        #np.power(10., a-ma.power(self.intersect, 1./2.)/ma.power(self.intersect, 1./2.)+np.log10(self.intersect)))

    def inverted(self):
        return SqrtLogTransform(intersect=self.intersect)




class SqrtLogLocator(Locator):
    """
    Determine the tick locations for symmetric log axes
    """

    def __init__(self, subs=(1.0,), sqrtthresh=None):
        """
        place ticks on the location= base**i*subs[j]
        """
        if sqrtthresh is not None:
            self._base = 10.
            self._sqrtthresh = sqrtthresh
        else:
            raise ValueError("linthresh must be provided.")
        if subs is None:
            self._subs = 'auto'
        else:
            self._subs = subs
        self.numticks = 15

    def set_params(self, subs=None, numticks=None):
        """Set parameters within this locator."""
        if numticks is not None:
            self.numticks = numticks
        if subs is not None:
            self._subs = subs

    def __call__(self):
        'Return the locations of the ticks'
        # Note, these are untransformed coordinates
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        b = self._base
        t = self._sqrtthresh

        logging.debug('vmin: %d, vmax: %d', vmin, vmax)

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        # The domain is divided into two sections, only some of
        # which may actually be present.
        #
        # 0== t ========>
        # aaa   bbbbbbbb
        #
        # b) will have ticks at integral log positions.  The
        # number of ticks needs to be reduced if there are more
        # than self.numticks of them.
        #
        # a) has major ticks corresponding to 0 and the two decades
        #    below t. Thus if t=10, major ticks will be at [0,0.1,1]
        #    Minor ticks at:
        #    [0.01, 0.05, 0.2, 0.4, 0.6, 0.8, 2, 3, 4, 5, 6, 7, 8, 9] * t/10
        #
        #    If t=4, major ticks will also be at [0,0.1,1]
        #
        # "simple" mode is when the range falls entirely within (0,
        # t) -- it should just display (0, vmax)

        has_a = has_b = False
        if vmin < t:
            has_a = True
            if vmax > t:
                has_b = True
        else:
            has_b = True

        def get_log_range(lo, hi):
            lo = np.floor(np.log(lo) / np.log(b))
            hi = np.ceil(np.log(hi) / np.log(b))
            return lo, hi

        # First, calculate all the ranges, so we can determine striding
        if has_b:
            if has_a:
                b_range = get_log_range(t, vmax + 1)
            else:
                b_range = get_log_range(vmin, vmax + 1)
        else:
            b_range = (0, 0)

        total_ticks = (b_range[1] - b_range[0])
        if has_b:
            total_ticks += 1
        stride = max(np.floor(float(total_ticks) / (self.numticks - 1)), 1)

        decades = []

        if has_a:
            decades.extend([0])
            decades.extend(b**(np.array([-2,-1]) + np.ceil(np.log(t)/np.log(b))))

        if has_b:
            decades.extend(b ** (np.arange(b_range[0], b_range[1], stride)))

        decades = np.unique(decades)
        
        # Add the subticks if requested

        if isinstance(self._subs, six.string_types):
            _first = 2.0 if self._subs == 'auto' else 1.0
            if len(decades) > 10 or b < 3:
                if self._subs == 'auto':
                    return np.array([])  # no minor or major ticks
                else:
                    subs = np.array([1.0])  # major ticks
            else:
                subs = np.arange(_first, b)
        else:
            subs = np.asarray(self._subs)

        if len(subs) > 1 or subs[0] != 1.0:
            ticklocs = []
            for decade in decades:
                if decade == 0: continue
                for sub in subs:
                    if sub * decade < vmax: ticklocs.append(sub * decade)
        else:
            ticklocs = decades

        logging.debug('b_range: %s', b_range)
        logging.debug('subs: %s', subs)
        logging.debug('decades: %s', decades)
        logging.debug('ticklocs: %s', ticklocs)

        return self.raise_if_exceeds(np.array(ticklocs))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'

        logging.debug('vmin: %s, vmax: %s', vmin, vmax)
        vmin, vmax = self.axis.get_data_interval()
        logging.debug('vmin: {}, vmax: {}', vmin, vmax)

        #pdb.set_trace()
        
        b = self._base
        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if True: #rcParams['axes.autolimit_mode'] == 'round_numbers':
            if vmin > np.log10(self._sqrtthresh)-2:
                #if not is_decade(abs(vmin), b):
                #    vmin = decade_down(vmin, b)
                vmin = 0.
            else:
                vmin = 0.
            if not is_decade(abs(vmax), b):
                if vmax < 0:
                    vmax = -decade_down(-vmax, b)
                else:
                    vmax = decade_up(vmax, b)

            if vmin == vmax:
                if vmin < np.log10(self._sqrtthresh):
                    vmin = 0.
                    vmax = np.log10(self._sqrtthresh)+1
                else:
                    vmin = decade_down(vmin, b)
                    vmax = decade_up(vmax, b)

        result = [vmin, vmax]
        logging.debug('adjusted: vmin: %s, vmax: %s', vmin, vmax)
        logging.debug('limits: %s', result)
        return result


        #'Try to choose the view limits intelligently'
        #b = self._base
        #if vmax < vmin:
        #    vmin, vmax = vmax, vmin
        #
        #if rcParams['axes.autolimit_mode'] == 'round_numbers':
        #    if not is_decade(abs(vmin), b):
        #        if vmin < 0:
        #            vmin = -decade_up(-vmin, b)
        #        else:
        #            vmin = decade_down(vmin, b)
        #    if not is_decade(abs(vmax), b):
        #        if vmax < 0:
        #            vmax = -decade_down(-vmax, b)
        #        else:
        #            vmax = decade_up(vmax, b)
        #
        #    if vmin == vmax:
        #        if vmin < 0:
        #            vmin = -decade_up(-vmin, b)
        #            vmax = -decade_down(-vmax, b)
        #        else:
        #            vmin = decade_down(vmin, b)
        #            vmax = decade_up(vmax, b)
        #
        #result = nonsingular(vmin, vmax)
        #return result

class SqrtLogScale(mscale.ScaleBase):
    """
    A standard logarithmic scale.  Care is taken so non-positive
    values are not plotted.

    For computational efficiency (to push as much as possible to Numpy
    C code in the common cases), this scale provides different
    transforms depending on the base of the logarithm:

       - base 10 (:class:`Log10Transform`)
       - base 2 (:class:`Log2Transform`)
       - base e (:class:`NaturalLogTransform`)
       - arbitrary base (:class:`LogTransform`)
    """
    name = 'sqrtlog'

    # compatibility shim

    def __init__(self, axis, **kwargs):
        """
        *basex*/*basey*:
           The base of the logarithm

        *nonposx*/*nonposy*: ['mask' | 'clip' ]
          non-positive values in *x* or *y* can be masked as
          invalid, or clipped to a very small positive number

        *subsx*/*subsy*:
           Where to place the subticks between each major tick.
           Should be a sequence of integers.  For example, in a log10
           scale: ``[2, 3, 4, 5, 6, 7, 8, 9]``

           will place 8 logarithmically spaced minor ticks between
           each major tick.
        """
        if axis.axis_name == 'x':
            intersect = kwargs.pop('intersectx', 1.)
            subs = kwargs.pop('subsx', None)
            nonpos = kwargs.pop('nonposx', 'mask')
        else:
            intersect = kwargs.pop('intersecty', 1.)
            subs = kwargs.pop('subsy', None)
            nonpos = kwargs.pop('nonposy', 'mask')

        if nonpos not in ['mask', 'clip']:
            raise ValueError("nonposx, nonposy kwarg must be 'mask' or 'clip'")

        self._transform = SqrtLogTransform(nonpos,intersect)

        self.intersect = intersect
        self.subs = subs

    def set_default_locators_and_formatters(self, axis):
        """
        Set the locators and formatters to specialized versions for
        log scaling.
        """
        #axis.set_major_locator(LogLocator(10.))
        #axis.set_major_formatter(LogFormatterSciNotation(10.))
        #axis.set_minor_locator(LogLocator(10., self.subs))
        #axis.set_minor_formatter(
        #    LogFormatterSciNotation(10.,
        #                            labelOnlyBase=(self.subs is not None)))

        axis.set_major_locator(SqrtLogLocator(sqrtthresh=self.intersect))
        axis.set_major_formatter(LogFormatterSciNotation(10.))
        axis.set_minor_locator(SqrtLogLocator(sqrtthresh=self.intersect, subs=self.subs))
        axis.set_minor_formatter(
            LogFormatterSciNotation(10.,
                                    labelOnlyBase=(self.subs is not None)))

    def get_transform(self):
        """
        Return a :class:`~matplotlib.transforms.Transform` instance
        appropriate for the given logarithm base.
        """
        return self._transform

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Limit the domain to positive values.
        """
        if not np.isfinite(minpos):
            minpos = 1e-300  # This value should rarely if ever
                             # end up with a visible effect.

        return (minpos if vmin < 0 else vmin,
                minpos if vmax < 0 else vmax)



def annotation_arrow( ax, xmin, xmax, y, text, linecolor='black', linewidth=1, fontsize=10, xycoords='data'):

#    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords=xycoords, textcoords=xycoords,
#            arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth,
#            'connectionstyle':'widthA=0.5,widthB=0.5'})
    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords=xycoords, textcoords=xycoords,
            arrowprops={'arrowstyle': '<->', 'color':linecolor, 'linewidth':linewidth}, zorder=-100)
    
    disp_xmin = ax.transData.transform((xmin, 0))[0]
    disp_xmax = ax.transData.transform((xmax, 0))[0]
    disp_xcenter = disp_xmin + (disp_xmax-disp_xmin)/2
    inv = ax.transData.inverted()
    xcenter = inv.transform((disp_xcenter, 0))[0]
    
    ax.annotate(text, xy=(xcenter,y), ha='center', va='center', fontsize=fontsize, 
                xycoords=xycoords, textcoords=xycoords,
                bbox=dict(fc='w', ec='w', pad=1), zorder=-99)
                
def annotate_xaxis(ax, intersect=None, arrows='top', arrow_margin=1.03):
    xmin = ax.get_xlim()[0]
    xmax = ax.get_xlim()[1]
    
    if intersect is not None:
        if arrows is not None or arrows == 'none':
            if arrows == 'top':
                annotation_arrow(ax,xmin,intersect,arrow_margin,'$\sqrt{t}$', xycoords=("data", "axes fraction"))
                annotation_arrow(ax,intersect,xmax,arrow_margin,'$\log{(t)}$', xycoords=("data", "axes fraction"))    
        ax.axvline(x=intersect, ls='--', color='k', lw=1)
                

                
# Now that the Scale class has been defined, it must be registered so
# that ``matplotlib`` can find it.
mscale.register_scale(SqrtLogScale)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    FORMAT = "[%(filename)s:%(funcName)s():%(lineno)s] %(message)s"
    logging.basicConfig(format=FORMAT, stream=sys.stderr, level=logging.DEBUG)

    f = plt.figure(figsize=(15,5))
    t = np.logspace(-1,3,50)
    s = np.log10(t)
    t2 = np.linspace(0, np.sqrt(1000), int((np.sqrt(1000)-np.sqrt(0))*10))**2
    s2 = np.sqrt(t2)
    intersect = 0.33

    plt.plot(t, s, '.-b', lw=2)
    plt.plot(t2, s2, '.-r', lw=2)
    #pdb.set_trace()

    ax = plt.gca()
    ax.set_xscale('sqrtlog', intersectx=intersect)
    ax.set_ylim([-2, np.ceil(np.sqrt(intersect))+1])
    #plt.gca().set_xlim([0, 2300])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('test')
    ax.grid(True)
    ax.grid(b=True, which='minor', ls='--')
    
    annotate_xaxis(ax=ax, intersect=intersect)
    
    #plt.show(block=False)