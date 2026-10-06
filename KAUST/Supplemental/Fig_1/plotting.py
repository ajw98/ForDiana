#!/usr/bin/env python3

# This program is used in conjunction with other programs that call for this program in order to make phase diagrams. This 
# specific version has been altered for making phase diagrams for Polysulfone, Polar Clean and Water

import math

import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np

"""Matplotlib Ternary plotting utility."""

## Constants ##

SQRT3OVER2 = math.sqrt(3) / 2.

## Default colormap, other options here: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
DEFAULT_COLOR_MAP = pyplot.get_cmap('jet')
#DEFAULT_COLOR_MAP = pyplot.get_cmap('gist_stern')

## Helpers ##
def unzip(l):
    return zip(*l)

def normalize(xs):
    """Normalize input list."""
    s = float(sum(xs))
    return [x / s for x in xs]

## Boundary ##

# old version of draw boundary
#def draw_boundary(scale=1.0, linewidth=2.0, color='black', ax=None):
#    # Plot boundary of 3-simplex.
#    if not ax:
#        ax = pyplot.subplot()
#    scale = float(scale)
#    # Note that the math.sqrt term is such to prevent noticable roundoff on the top corner point.
#    ax.plot([0, scale, scale / 2, 0], [0, 0, math.sqrt(scale * scale * 3.) / 2, 0], color, linewidth=linewidth)
#    ax.set_ylim((-0.05 * scale, .90 * scale))
#    ax.set_xlim((-0.05 * scale, 1.05 * scale))
#    return ax

# new version edited by DRT on 01/29/15
def draw_boundary(scale=1.0, linewidth=1.5, color='black', ax=None):
    # Plot boundary of 3-simplex.
    if not ax:
        ax = pyplot.subplot()
    scale = float(scale)
    # Note that the math.sqrt term is such to prevent noticable roundoff on the top corner point.
    ax.plot([0, scale, scale / 2, 0], [0, 0, math.sqrt(scale * scale * 3.) / 2, 0], color, linewidth=linewidth)
    # set up axes
    ax.set_xlim((-0.1 * scale, 1.1 * scale))
    ax.set_ylim((-0.1 * scale, 1.0 * scale))
    # set up axis labels
    #ax.text(-0.05, -0.03, "$\phi_{1}$", fontsize=20)
    #ax.text( 1.02, -0.03, "$\phi_{2}$", fontsize=20)
    #ax.text( 0.48, math.sqrt(3.)/2.+0.03, "$\phi_{3}$", fontsize=20)

#    # set up lines on plot
    #dline = 0.1
    #n = int((1./dline)-1)
    #lines = np.linspace(0, 1, n+2)[1:-1]
    #phi1_L = np.vstack([lines, np.zeros(n), 1-lines]).T
    #phi1_R = np.vstack([lines, 1-lines, np.zeros(n)]).T
    #phi2_L = np.vstack([np.zeros(n), lines, 1-lines]).T
    #phi2_R = np.vstack([1-lines, lines, np.zeros(n)]).T
    #phi3_L = np.vstack([np.zeros(n), 1-lines, lines]).T
    #phi3_R = np.vstack([1-lines, np.zeros(n), lines]).T

    #[phi1_Lx, phi1_Ly] = project(phi1_L)
    #[phi1_Rx, phi1_Ry] = project(phi1_R)
    #[phi2_Lx, phi2_Ly] = project(phi2_L)
    #[phi2_Rx, phi2_Ry] = project(phi2_R)
    #[phi3_Lx, phi3_Ly] = project(phi3_L)
    #[phi3_Rx, phi3_Ry] = project(phi3_R)

    #ax.plot([phi1_Lx, phi1_Rx], [phi1_Ly, phi1_Ry], 'k-', linewidth=0.5)
    #ax.plot([phi2_Lx, phi2_Rx], [phi2_Ly, phi2_Ry], 'k-', linewidth=0.5)
    #ax.plot([phi3_Lx, phi3_Rx], [phi3_Ly, phi3_Ry], 'k-', linewidth=0.5)
#    ax.set_facecolor("white")
    ax.axes.get_xaxis().set_visible(True)
    #ax.axes.get_yaxis().set_visible(False)
    ax.text(-0.07, -0.07, 'Cellulose Acetate', fontsize=10)
    ax.text( .9, -0.07, 'Water', fontsize=10)
    ax.text( 0.27, math.sqrt(3.)/2.+0.03, 'Glacial Acetic Acid', fontsize=10)
    #ax.text(0.22,0.15,'$2\phi$',fontsize=12)
    #ax.text(0.5,0.4,'$\mu\phi$', fontsize=12)

    # Text for axis
    # ax.text(0.1,-0.05,'0.10')



    return ax



## Curve Plotting ##
def project_point(p):
    """Maps (x,y,z) coordinates to planar-simplex."""
    a = p[0]
    b = p[1]
    c = p[2]
    x = 0.5 * (2 * b + c)
    y = SQRT3OVER2 * c
    return (x, y)

def project(s):
    """Maps (x,y,z) coordinates to planar-simplex."""
    # Is s an appropriate sequence or just a single point?
    try:
        return unzip(map(project_point, s))
    except TypeError:
        return project_point(s)
    except IndexError: # for numpy arrays
        return project_point(s)

def plot(t, color=None, linewidth=1.0, ax=None):
    """Plots trajectory points where each point satisfies x + y + z = 1. First argument is a list or numpy array of tuples of length 3."""
    if not ax:
        ax = pyplot.subplot()
    xs, ys = project(t)
    if color:
        ax.plot(xs, ys, c=color, linewidth=linewidth)
    else:
        ax.plot(xs, ys, linewidth=linewidth)
    return ax


## Heatmaps##

def simplex_points(steps=100, boundary=True):
    """Systematically iterate through a lattice of points on the 2 dimensional simplex."""
    steps = steps - 1
    start = 0
    if not boundary:
       start = 1
    for x1 in range(start, steps + (1-start)):
        for x2 in range(start, steps + (1-start) - x1):
            x3 = steps - x1 - x2
            yield (x1, x2, x3)

def colormapper(x, a=0, b=1, cmap=None):
    """Maps color values to [0,1] and obtains rgba from the given color map for triangle coloring."""
    if b - a == 0:
        rgba = cmap(0)
    else:
        rgba = cmap((x - a) / float(b - a))
    hex_ = matplotlib.colors.rgb2hex(rgba)
    return hex_

def triangle_coordinates(i, j, alt=False):
    """Returns the ordered coordinates of the triangle vertices for i + j + k = N. Alt refers to the averaged triangles; the ordinary triangles are those with base parallel to the axis on the lower end (rather than the upper end)"""
    # N = i + j + k
    if not alt:
        return [(i/2. + j, i * SQRT3OVER2), (i/2. + j + 1, i * SQRT3OVER2), (i/2. + j + 0.5, (i + 1) * SQRT3OVER2)]
    else:
        # Alt refers to the inner triangles not covered by the default case
        return [(i/2. + j + 1, i * SQRT3OVER2), (i/2. + j + 1.5, (i + 1) * SQRT3OVER2), (i/2. + j + 0.5, (i + 1) * SQRT3OVER2)]

def heatmap(d, steps, cmap_name=None, boundary=True, ax=None, scientific=False):
    """Plots values in the dictionary d as a heatmap. d is a dictionary of (i,j) --> c pairs where N = steps = i + j + k."""
    if not ax:
        ax = pyplot.subplot()
    if not cmap_name:
        cmap = DEFAULT_COLOR_MAP
    else:
        cmap = pyplot.get_cmap(cmap_name)
    a = min(d.values())
    b = max(d.values())
    # Color data triangles.
    for k, v in d.items():
        i, j = k
        vertices = triangle_coordinates(i,j)
        x,y = unzip(vertices)
        color = colormapper(d[i,j],a,b,cmap=cmap)
        ax.fill(x, y, facecolor=color, edgecolor=color)
    # Color smoothing triangles.
    offset = 0
    if not boundary:
        offset = 1
    for i in range(offset, steps+1-offset):
        for j in range(offset, steps -i -offset):
            try:
                alt_color = (d[i,j] + d[i, j + 1] + d[i + 1, j])/3.
                color = colormapper(alt_color, a, b, cmap=cmap)
                vertices = triangle_coordinates(i,j, alt=True)
                x,y = unzip(vertices)
                pyplot.fill(x, y, facecolor=color, edgecolor=color)
            except KeyError:
                # Allow for some portions to have no color, such as the boundary
                pass
    # Colorbar hack
    # http://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
    sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=pyplot.Normalize(vmin=a, vmax=b))
    # Fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cb = pyplot.colorbar(sm, ax=ax, format='%.3f')
    if scientific:
        cb.formatter = matplotlib.ticker.ScalarFormatter()
        cb.formatter.set_powerlimits((0, 0))
        cb.update_ticks()
    return ax

## Convenience Functions ##
    
def plot_heatmap(func, steps=40, boundary=True, cmap_name=None, ax=None):
    """Computes func on heatmap coordinates and plots heatmap. In other words, computes the function on sample points of the simplex (normalized points) and creates a heatmap from the values."""
    d = dict()
    for x1, x2, x3 in simplex_points(steps=steps, boundary=boundary):
        d[(x1, x2)] = func(normalize([x1, x2, x3]))
    heatmap(d, steps, cmap_name=cmap_name, ax=ax)
    
def plot_multiple(trajectories, linewidth=2.0, ax=None):
    """Plots multiple trajectories and the boundary."""
    if not ax:
        ax = pyplot.subplot()
    for t in trajectories:
        plot(t, linewidth=linewidth, ax=ax)
    draw_boundary(ax=ax)
    return ax
