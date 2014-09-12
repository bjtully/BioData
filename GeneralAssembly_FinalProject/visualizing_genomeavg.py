#!/usr/bin/python
"""This script takes the output of multiGenomeARSC_V2.py - a tab-delimited text file with the ARSC calculations -
and plots them using the Bokeh visualization module.

USAGE: python visualizing_genomeavg.py <precalculated_avgARSC.txt> <ELEMENT OF CHOICE: nitrogen, sulfur, carbon, oxygen>"""

from __future__ import division
import pylab as pl
import sys

import numpy as np
from six.moves import zip
from collections import OrderedDict
from bokeh.plotting import *
from bokeh.objects import HoverTool
from bokeh.objects import Range1d
import pandas
from patsy import dmatrices
from sklearn.linear_model import LinearRegression

#Read precalculated_avgARSC.txt into a data format
try:
	data = pandas.read_csv("%s" % str(sys.argv[1]), sep="\t")
except IndexError:
	print "Please provide correct input.\nUsage: python visualizing_genomeavg.py <precalculated_avgARSC.txt> <ELEMENT OF CHOICE: nitrogen, sulfur, carbon, oxygen>"
	exit()
except IOError:
	print "Please provide correct input.\nUsage: python visualizing_genomeavg.py <precalculated_avgARSC.txt> <ELEMENT OF CHOICE: nitrogen, sulfur, carbon, oxygen>"
	exit()
#set output file name
try:
	output_file("%s.html" % str(sys.argv[2]).lower())
except IndexError:
	print "Please provide correct input.\nUsage: python visualizing_genomeavg.py <precalculated_avgARSC.txt> <ELEMENT OF CHOICE: nitrogen, sulfur, carbon, oxygen>"
	exit()

#set empty list sets for storing data
accession = []
species_name = []
X = []
n_y = []
c_y = []
s_y = []
o_y = []

#parse through the tab delimited text file - order:
#		0 				1 			  2 		3     		4 			  5 		   6
#Accession ID | Full organism Name | %GC | Avg N ARSC | Avg C ARSC | Avg S ARSC | Avg O ARSC
for z in open(str(sys.argv[1]), "r"):
    if z[:9] != "Accession":
        z = z.rstrip()
        a = z.split("\t")
        accession.append(a[0])
        species_name.append(a[1])
        X.append(float(a[2]))
        n_y.append(float(a[3]))
        c_y.append(float(a[4]))
        s_y.append(float(a[5]))
        o_y.append(float(a[6]))

#following the Bokeh tutorial about Hover Tools
#set tools to be used by file graph
TOOLS = "pan,wheel_zoom,box_zoom,reset,previewsave,hover"
#the inds variable is needed for the text parameter below. it stores an "index" value counting up for each element
#in the list
N = len(X)
inds = [str(i) for i in np.arange(N)]

#following the Bokeh tutorial about Hover Tools
source = ColumnDataSource(
    data=dict(
        x=X,
        y=c_y,
        accession=accession,
        species_name=species_name,
    )
)
#information needed to plot regression lines
#points plotted along regression line
n = range(0,110,10)
#empty list to store y values
regression_y = []

#Nitrogen. With more time some of this could be condensed into a function
if str(sys.argv[2]).lower() == "nitrogen":
	#setting the max values from the X values (X aka %GC) and Y values (Avg N ARSC)
	min_x = min(X)
	max_x = max(X)
	min_y = min(n_y)
	max_y = max(n_y)
	#setting range as detailed in Bokeh tutorial
	xr = Range1d(start=min_x-2, end=max_x+2)
	yr = Range1d(start=min_x-0.5, end=max_x+0.5)
	#hold required for showing both the scatter points and line
	hold()
	figure(
		title="Average Nitrogen Atom Per Residue Side Chain",
		)
	#Bokeh scatter item
	#X = x-values
	#n_y = y-values
	#source set above. tools set above
	#high alpha for solid circles. with pt = 5
	scatter(X, n_y, size = 5, source=source, 
	        color = "gray", alpha = 0.9, tools=TOOLS)
	#Bokeh Hover text item
	#inds contains an "index" value without any meaning, set alpha to 0 so they are 100% transparent
	#standard text. aligned center. on the baseline. horizontal text.
	text(X, n_y, text=inds, alpha=0, 
	     text_font_size="5pt", text_baseline="middle", 
	     text_align="center", angle=0)
	#plotting the regression libe
	#calculating the y-values for this particular line based on a linear regression model
	#the linear regression model accepts the %GC and determines Avg ARSC
	y, X = dmatrices('AvgN_ARSC ~ GC - 1', data=data, return_type="dataframe")
	
	model = LinearRegression()
	model = model.fit(X,y)
	coef = model.coef_
	intercept = model.intercept_
	for i in n:
		regression_y.append(coef*i + intercept)
	#Bokeh line item
	line(n, regression_y, color="black", line_width=2)
	#display R^2 value
	print "The R^2 value is:\t"+str(model.score(X,y))
	#following Bokeh tutorial. creates the hover tool
	hover = [t for t in curplot().tools if isinstance(t, HoverTool)][0]
	#information displayed in the Hover item
	hover.tooltips = OrderedDict([
	    ("GC content", "@x"),
	    ("Avg. N ARSC", "@y"),
	    ("Species Name", "@species_name"),
	    ("Accession ID", "@accession")
	])
	show()
#Identical as Nitrogen but for Carbon
if str(sys.argv[2]).lower() == "carbon":
	min_x = min(X)
	max_x = max(X)
	min_y = min(c_y)
	max_y = max(c_y)

	xr = Range1d(start=min_x-2, end=max_x+2)
	yr = Range1d(start=min_x-0.5, end=max_x+0.5)

	hold()
	figure(
		title="Average Carbon Atom Per Residue Side Chain",
		)
	scatter(X, c_y, size = 5, source=source, 
	        color = "green", alpha = 0.9, tools=TOOLS)
	text(X, c_y, text=inds, alpha=0, 
	     text_font_size="5pt", text_baseline="middle", 
	     text_align="center", angle=0)
	y, X = dmatrices('AvgC_ARSC ~ GC - 1', data=data, return_type="dataframe")
	
	model = LinearRegression()
	model = model.fit(X,y)
	coef = model.coef_
	intercept = model.intercept_
	for i in n:
		regression_y.append(coef*i + intercept)

	line(n, regression_y, color="purple", line_width=2)
	print "The R^2 value is:\t"+str(model.score(X,y))
	hover = [t for t in curplot().tools if isinstance(t, HoverTool)][0]
	hover.tooltips = OrderedDict([
	    ("GC content", "@x"),
	    ("Avg. C ARSC", "@y"),
	    ("Species Name", "@species_name"),
	    ("Accession ID", "@accession")
	])
	show()
#Identical as Nitrogen but for Sulfur
if str(sys.argv[2]).lower() == "sulfur":
	min_x = min(X)
	max_x = max(X)
	min_y = min(s_y)
	max_y = max(s_y)

	xr = Range1d(start=min_x-2, end=max_x+2)
	yr = Range1d(start=min_x-0.5, end=max_x+0.5)

	hold()
	figure(
		title="Average Sulfur Atom Per Residue Side Chain",
		)
	scatter(X, s_y, size = 5, source=source, 
	        color = "blue", alpha = 0.9, tools=TOOLS)
	text(X, s_y, text=inds, alpha=0, 
	     text_font_size="5pt", text_baseline="middle", 
	     text_align="center", angle=0)
	y, X = dmatrices('AvgS_ARSC ~ GC - 1', data=data, return_type="dataframe")
	
	model = LinearRegression()
	model = model.fit(X,y)
	coef = model.coef_
	intercept = model.intercept_
	for i in n:
		regression_y.append(coef*i + intercept)
	
	line(n, regression_y, color="gray", line_width=2)
	print "The R^2 value is:\t"+str(model.score(X,y))
	hover = [t for t in curplot().tools if isinstance(t, HoverTool)][0]
	hover.tooltips = OrderedDict([
	    ("GC content", "@x"),
	    ("Avg. S ARSC", "@y"),
	    ("Species Name", "@species_name"),
	    ("Accession ID", "@accession")
	])
	show()
#Identical as Nitrogen but for Oxygen
if str(sys.argv[2]).lower() == "oxygen":
	min_x = min(X)
	max_x = max(X)
	min_y = min(o_y)
	max_y = max(o_y)

	xr = Range1d(start=min_x-2, end=max_x+2)
	yr = Range1d(start=min_x-0.5, end=max_x+0.5)

	hold()
	figure(
		title="Average Oxygen Atom Per Residue Side Chain",
		)
	scatter(X, o_y, size = 5, source=source, 
	        color = "orange", alpha = 0.9, tools=TOOLS)
	text(X, o_y, text=inds, alpha=0, 
	     text_font_size="5pt", text_baseline="middle", 
	     text_align="center", angle=0)
	y, X = dmatrices('AvgO_ARSC ~ GC - 1', data=data, return_type="dataframe")
	
	model = LinearRegression()
	model = model.fit(X,y)
	coef = model.coef_
	intercept = model.intercept_
	for i in n:
		regression_y.append(coef*i + intercept)
	
	line(n, regression_y, color="blue", line_width=2)
	print "The R^2 value is:\t"+str(model.score(X,y))
	hover = [t for t in curplot().tools if isinstance(t, HoverTool)][0]
	hover.tooltips = OrderedDict([
	    ("GC content", "@x"),
	    ("Avg. O ARSC", "@y"),
	    ("Species Name", "@species_name"),
	    ("Accession ID", "@accession")
	])
	show()