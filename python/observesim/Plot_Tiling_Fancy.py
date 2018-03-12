import iris
import pandas as pd
import numpy as np
import holoviews as hv
import geoviews as gv
import geoviews.feature as gf 
from cartopy import crs
#hv.notebook_extension()

test = pd.read_csv('/Users/jsobeck/observesim/data/tiling/test.csv')
test = hv.Dataset(test)
test.data.head()

output size=200
(features.relabel(group='Mollweide')(plot=dict(projection=crs.Mollweide()))

%%opts Feature.Lines [projection=ccrs.Mollweide()] (facecolor='none' edgecolor='gray')

dataset = gv.Dataset(test, kdims=['lon', 'lat'], crs=crs.Mollweide())
dataset.to(gv.Image, ['lon', 'lat'])
