# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:31:53 2022

This module includes tools for analysing strings.

@author: Ketil
"""
from distutils.util import strtobool

def strornum(string):
    """
    Converts a string to integer or float. If the string is neither, return the string.

    Parameters
    ----------
    string : str
             String containing some numerical value. If the string is not either an int or a float, return the pure string.

    Returns
    -------
    int, float or str
             Either an integer, float or a string.
    """
    try:
        return int(string)
    except ValueError:
        pass
    try:
        return float(string)
    except ValueError:
        return str(string)


def strtolist(string, dtype=str):
    """
    Construct a list out of a comma separated string. Blank spaces are removed between elements. 
    Set the type of input to be str, int, float etc.

    Parameters
    ----------
    string : str
             String to be parsed into a list.
    dtype  : dtype, optional
             Type of each element of list, by default str.

    Returns
    -------
    list
        List of elements separated by comma in original string.
    """
    return list(map(dtype, string.strip('[]').replace(' ','').replace("\'","").split(',')))


def strtoBool(string):
    """
    Returns a Boolean based on the string input.
    If string = 'True' then it returns True.
    If string = 'False' then it returns False.

    Parameters
    ----------
    string : str
             String containing either 'True' or 'False'

    Returns
    -------
    bool
        True or False
    """
    return bool(strtobool(string))


def any_digits(string):
    """
    Returns True if string contains any digits. Returns False if string contains no digits.

    Parameters
    ----------
    string : str
             Some string.

    Returns
    -------
    bool
        True if string contains digits, False if string contains no digits.
    """
    return any(i.isdigit() for i in string)