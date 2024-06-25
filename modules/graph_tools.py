import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def name_plot_corner(ax, text, weight = "bold", size = "x-large", x = 0.02, 
                     y = 0.92, **txtkwargs):
    """Add a name in the corner of a plot.

    Args:
        ax (matplotlib Axes object): ax on which to plot the text 
        text (str): text to plot
        weight (str, optional): weight of the text. Defaults to "bold".
        size (str, optional): size of the text. Defaults to "x-large".
        x (float, optional): position of the text. Defaults to 0.02.
        y (float, optional): position of the text. Defaults to 0.92.
    Notes:
        positions of the text are in axis coordinates (not in data coordinates).
    """
    ax.text(x = x, y = y, s = text, weight = weight, 
             size = size, transform = ax.transAxes, **txtkwargs)

def mix_color(c1,c2,c3,u1,u2,u3, geom = True):
    """ Mixes 3 RGB color according to weights.

    Args:
        c1, c2, c3 (arrays): size-3 arrays containing the RGB colors
        u1, u2, u3 (floats): non-negative weights, at least one must be > 0.
        geom (bool, optional): Geometric average. If False, the color are mixed 
           with arithmetic average. Defaults to True.

    Returns:
        array: mixed color in RGB space.
    """
    if geom:# Geometric mean
        return np.power(c1**u1 * c2**u2 * c3**u3, 1/(u1+u2+u3))
    return (c1*u1 + c2*u2 + c3*u3)/(u1+u2+u3)

def plot_grad_line(ax, x0, y0, x1, y1, n_ticks = 5,                   
                   lsty = dict(ls = "-", color = "k"), 
                   revert_marker = False, ticksize = 3.0):
    """Plot a graduated line on an axis

    Args:
        ax (matplotlib Axes object): ax on which to plot the line.
        x0, y0, x1, y1 (floats): coordinates of the start and end of the line.
        n_ticks (int, optional): number of ticks. Defaults to 5.
        lsty (dict, optional): line style, kwargs accepted by plot function. 
            Defaults to dict(ls = "-", color = "k").
        revert_marker (bool, optional): allows to change the side of the line on
            which to draw the ticks. Defaults to False.
        ticksize (float, optional): Tick size. Defaults to 3.0
    """
    X = np.linspace(x0, x1, n_ticks)
    Y = np.linspace(y0, y1, n_ticks)
    m = mpl.markers.MarkerStyle(3)
    angle = np.pi * float(revert_marker) + np.arctan((y1-y0) / (x1-x0))
    m._transform.rotate(angle)
    ax.plot(X, Y, marker = m, **lsty, markersize = ticksize)


def color_by_distance(xx, yy, xtri, ytri, col1, col2, col3, geom):
    # Mix color according to the position of a point (xx, yy) with respect to 
    # the vertices of a triangle of coordinates (-xtri, 0), (xtri, 0), (0, ytri).
    d1 = 1-np.sqrt((xx - xtri)**2 + yy**2)
    d2 = 1-np.sqrt(xx**2 + (yy - ytri)**2)
    d3 = 1-np.sqrt((xx+xtri)**2 + yy**2)
    return mix_color(col1, col2, col3, d1, d2, d3, geom = geom)

def color_scale_triangle(ax, col1, col2, col3, geom = True):
    """ Plot a color scale with 3 colors with the shape of a triangle.

    Args:
        ax (matplotlib Axes object): ax on which to draw the triangle.
        col1, col2, col3 (arrays): size-3 arrays containing the RGB colors
        geom (bool, optional): Geometric average. If False, the color are mixed 
           with arithmetic average. Defaults to True.
    """
    ax.axis("off")
    xtri, ytri = 0.5, np.sqrt(0.75) # coordinates triangle

    x = np.linspace(-xtri-0.25, xtri + 0.25, 300)
    y = np.linspace(-0.25, ytri + 0.25, 300)
    captionmat = np.ones((len(y), len(x), 4)) # RGB - alpha
    captionmat[:,:,3] = 0.0
    for j in range(len(y)):
        yy = y[j]
        if yy>0 and yy <= 1:
            for i in range(len(x)):
                xx = x[i]
                if yy <= ytri - np.abs(xx)*ytri/xtri:
                    captionmat[j,i,:3] = color_by_distance(xx, yy, xtri, ytri, col1, col2, col3, geom)
                    captionmat[j,i,3] = 1.0 # opaque

    ax.pcolormesh(x, y, captionmat, rasterized = True)
    plot_grad_line(ax, -xtri, 0, xtri, 0, 6)
    plot_grad_line(ax, -xtri, 0, 0, ytri, 6, revert_marker=True)
    plot_grad_line(ax, xtri, 0, 0, ytri, 6, revert_marker=True)