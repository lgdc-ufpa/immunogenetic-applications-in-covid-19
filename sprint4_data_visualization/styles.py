import seaborn as sns
import matplotlib.pyplot as plt
from typing import Optional


################
# Seaborn - Grid
################

def sns_style_grid_darkgrid():
    sns.set_style('darkgrid')


def sns_style_grid_white():
    sns.set_style('white')


def sns_style_grid_dark():
    sns.set_style('dark')


def sns_style_grid_whitegrid():
    sns.set_style('whitegrid')


def sns_style_grid_ticks():
    sns.set_style('ticks')


#################
# Seaborn - Spine
#################


def sns_style_spine_remove(
        left: Optional[bool] = False,
        **kwargs) -> None:
    """despine function

    Arguments
    ---------

    left : left only bottom spine

    kwargs : keyword arguments

    -> return : None
    """
    sns.despine(left=left, **kwargs)


###################
# Seaborn - context
###################


def sns_context(context=None, font_scale=1, rc=None) -> None:
    sns.set_context(context=context, font_scale=font_scale, rc=rc)


def sns_context_context_poster() -> None:
    sns.set_context(context='poster')


def sns_context_context_poster_font_scale_2() -> None:
    sns.set_context(context='poster', font_scale=2)
#####################
# Matplotlib - Figure
#####################


def plt_figure(**kwargs) -> None:
    plt.figure(**kwargs)


def plt_figure_size(x, y) -> None:
    plt.figure(figsize=(x, y))
