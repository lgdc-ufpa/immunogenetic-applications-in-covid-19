import seaborn as sns
import matplotlib.pyplot as plt
from typing import Optional, List, Dict, Union
from abc import ABC, abstractmethod


####################
# Distribution plots
####################


def distplot(
	df, 
	column: str, 
	*args, 
	**kwargs) -> sns.displot:
	"""Distplot function

	Arguments
	---------

	df : pandas dataframe

	columns: column of df

	args: arguments

	kwargs: keyword arguments

	-> return : sns.distplot
	"""

	return sns.displot(
		df[column], 
		*args, 
		**kwargs)


def jointplot(
	df, 
	x: str, 
	y: str, 
	*args, **kwargs) -> sns.jointplot:
	"""jointplot function

	Arguments
	---------

	df : pandas dataframe

	x : column in df

	y : column in df

	args: arguments

	kwargs: keyword arguments

	-> return : sns.jointplot
	"""

	return sns.jointplot(x=x, y=y, data=df, *args, **kwargs)
	

def pairplot(
	df, 
	hue: Optional[str] = None, 
	palette: Optional[str] = None, 
	*args, **kwargs) -> sns.pairplot:
	"""pairplot function

	Arguments
	---------
	
	df : pandas dataframe

	hue : column in df

	palette : color scheme

	-> return : sns.pairplot
	"""

	return sns.pairplot(
		df,
		hue=hue,
		palette=palette,
		*args,
		**kwargs)


def rugplot(
	df, 
	column: str, 
	*args, **kwargs) -> sns.rugplot:
	"""rugplot function

	Arguments
	---------

	df : pandas dataframe

	column : column in df

	args : arguments

	kwargs : keyword arguments

	-> return : sns.rugplot

	"""
	return sns.rugplot(
		df[column], 
		*args, 
		**kwargs)


def kdeplot(
	df, 
	column: str, 
	*args, **kwargs) -> sns.kdeplot:
	"""kdeplot function

	Arguments
	---------

	df : pandas dataframe

	column : column in df

	args : arguments

	kwargs : keyword arguments

	-> return : sns.kdeplot
	"""

	return sns.kdeplot(
		df[column],
		*args,
		**kwargs)


###################
# Categorical plots
###################


def barplot(
	df, 
	x: Optional[str] = None, 
	y: Optional[str] = None,	
	*args, **kwargs) -> sns.barplot:
	"""barplot function

	Arguments
	---------

	df : pandas dataframe

	x : column in dataframe

	y : column in dataframe

	args : arguments

	kwargs : keyword arguments

	-> return sns.barplot

	"""
	return sns.barplot(
		x=x,
		y=y,
		data=df,
		*args,
		**kwargs)	


def countplot(
	df,
	x: str,
	*args, **kwargs) -> sns.countplot:
	"""countplot function

	Arguments
	---------

	df : pandas dataframe

	x : column in df
	
	args : arguments

	kwargs : keyword arguments

	-> return : sns.countplot
	"""

	return sns.countplot(
		data=df,
		x=x,
		*args,
		**kwargs)


def boxplot(
	df,
	x: Optional[str] = None,
	y: Optional[str] = None,
	palette: Optional[str] = None,
	hue: Optional[str] = None,
	orient: Optional[str] = None,
	*args, **kwargs) -> sns.boxplot:
	"""boxplot function

	Arguments
	---------

	df : pandas dataframe

	x : column in df

	y : column in df

	palette : color scheme

	hue : column in df

	orient: v (vertical) or h (horizontal)

	args : arguments

	kwargs : keyword arguments

	-> return : sns.boxplot
	"""
	
	return sns.boxplot(
		x=x,
		y=y,
		data=df,
		palette=palette,
		hue=hue,
		orient=orient,
		*args,
		**kwargs)


def violinplot(
	df,
	x: str,
	y: str,
	hue: Optional[str] = None,
	split: Optional[bool] = False,
	*args,
	**kwargs) -> sns.violinplot:
	"""violinplot function
	
	Arguments
	---------

	df : pandas dataframe

	x : column in df

	y : column in df

	hue : column in df

	args : arguments

	kwargs : keyword arguments

	-> return : sns.violinplot
	"""
	return sns.violinplot(
		data=df,
		x=x,
		y=y,
		hue=hue,
		split=split,
		*args,
		**kwargs)


def stripplot(
	df,
	x: str,
	y: str,
	jitter: Optional[bool] = False,
	hue: Optional[str] = None,
	dodge: Optional[bool] = False,
	*args, **kwargs) -> sns.stripplot:
	"""stripplot function
	
	Arguments
	---------

	df : pandas dataframe

	x : column in df

	y : column in df

	jitter : to sparse points, or not

	hue : column in df

	dodge : to sparse points (util when use hue)

	args : arguments

	kwargs : keyword arguments

	-> return : sns.stripplot
	"""
	return sns.stripplot(
		data=df,
		x=x,
		y=y,
		jitter=jitter,
		hue=hue,
		dodge=dodge,
		*args,
		**kwargs)


def swarmplot(**kwargs) -> sns.swarmplot:
	"""swarmplot function

	Arguments
	---------

	data : pandas dataframe

	x : column in df

	y : column in df

	hue : column in df	

	dodge : to sparse points (util when use hue)

	color : color of points

	kwargs : keyword arguments

	-> return : sns.swarmplot
	"""

	return sns.swarmplot(**kwargs)


def catplot(
	df,
	x: str,
	y: str,
	kind: Optional[str] = None,
	**kwargs) -> sns.catplot:
	"""factorplot function
	
	Arguments
	---------

	df : pandas dataframe

	x : column in df

	y : column in df

	kind : kind of plot. Egg: 'bar'

	-> return : sns.catplot
	"""

	return sns.catplot(
		data=df,
		x=x,
		y=y,
		kind=kind,
		**kwargs)


#################
# Matricial Plots
#################


def heatmap(
	data,
	cmap: Optional[str] = None,
	annot: Optional[str] = None,
	linecolor: Optional[str] = None,
	linewidths: Optional[float] = None,
	**kwargs) -> sns.heatmap:
	"""heatmap function

	Arguments
	---------

	data : dataframe

	cmap : color scheme

	annot : number annotations

	linecolor : color of the line. Egg: 'white', '#FFFFFF'

	linewidths : width of the lines

	kwargs : keyword arguments

	-> return : sns.heatmap
	"""

	return sns.heatmap(
		data=data,
		cmap=cmap,
		annot=annot,
		linecolor=linecolor,
		linewidths=linewidths,
		**kwargs)


def clustermap(
	data,
	standard_scale: Optional[float] = None,
	**kwargs) -> sns.clustermap:
	"""clustermap function

	Arguments
	---------
	
	data : pandas dataframe

	standard_scale : scale of the colors of the clustermap

	kwargs : keyword arguments

	-> return : sns.clustermap

	"""
	return sns.clustermap(
		data=data,
		standard_scale=standard_scale,
		**kwargs)


##################
# Regression plots
##################


def lmplot(
	data,
	x: str,
	y: str,
	hue: Optional[str] = None,
	palette: Optional[str] = None,	
	markers: Optional[List[str]] = "o",
	scatter_kws: Optional[Dict] = None,
	col: Optional[str] = None,
	row: Optional[str] = None,
	aspect: Optional[float] = 1,
	height: Optional[float] = None,
	**kwargs) -> sns.lmplot:
	"""linear model plot function

	Arguments
	---------

	data : pandas dataframe

	x : column in pandas df

	y : column in pandas df

	hue : column in df

	palette : color scheme

	markers : point markers. Egg: ['o'], ['o', 'v']

	scatter_kws : keyword scatter arguments. Egg: {'s':100}

	col : column in df
	Used to segregate/split in columns the plot by selected column items
	Utils when used with hue. Egg: hue='sex', col='sex'

	row : column in df
	Used to segregate/split in rows the plot by selected column items
	Utisl when used with hue

	aspect : alter plot horizontal layout

	height : alter plot vertical layout

	kwargs : keyword arguments

	-> return : sns.lmplot
	"""
	return sns.lmplot(
		data=data,
		x=x,
		y=y,
		hue=hue,
		palette=palette,
		markers=markers,
		scatter_kws=scatter_kws,
		col=col,
		row=row,
		aspect=aspect,
		height=height,
		**kwargs)


###############
# Pairgrid plot
###############


class IPairGrid(ABC):
	def __init__(
		self,
		data,
		plot_type: Union[plt.scatter, plt.hist, sns.kdeplot],
		**kwargs):

		self.pairgrid = sns.PairGrid(data, **kwargs)
		self.plot_type = plot_type

	@abstractmethod
	def plot(self):
		pass	
		

class PairGridMap(IPairGrid):
	def __init__(self, data, plot_type=plt.scatter, **kwargs):
		IPairGrid.__init__(self, data, plot_type, **kwargs)

	def plot(self):
		self.pairgrid.map(self.plot_type)


class PairGridMapDiag(IPairGrid):
	def __init__(self, data, plot_type=plt.hist, **kwargs):
		IPairGrid.__init__(self, data, plot_type, **kwargs)

	def plot(self):
		return self.pairgrid.map_diag(self.plot_type)


class PairGridMapUpper(IPairGrid):
	def __init__(self, data, plot_type=plt.scatter, **kwargs):
		IPairGrid.__init__(self, data, plot_type, **kwargs)

	def plot(self):
		return self.pairgrid.map_upper(self.plot_type)


class PairGridMapLower(IPairGrid):
	def __init__(self, data, plot_type=sns.kdeplot, **kwargs):
		IPairGrid.__init__(self, data, plot_type, **kwargs)

	def plot(self):
		return self.pairgrid.map_lower(self.plot_type)


###############
# Facetgrid plot
###############


class IFacetGrid(ABC):
	def __init__(
		self,
		data,
		col: str,
		row: str,
		plot_type: Union[plt.scatter, plt.hist, sns.kdeplot],		
		**kwargs):

		self.facetgrid = sns.FacetGrid(data=data, col=col, row=row, **kwargs)
		self.plot_type = plot_type

	@abstractmethod
	def plot(self):
		pass


class FacetGridMap(IFacetGrid):
	def __init__(self, data, col, row, plot_type=plt.hist, **kwargs):
		IFacetGrid.__init__(self, data, col, row, plot_type, **kwargs)

	def plot(self, y: str):
		return self.facetgrid.map(self.plot_type, y)
