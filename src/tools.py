import matplotlib.pylab as pl

def remove_top_right_axis(ax=None):
    if ax is None:
        ax = pl.gca()

    # Remove top and right axis of parent ax
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().tick_left()
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()

