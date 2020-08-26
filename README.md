# [genice-twist](https://github.com/vitroid/genice-twist/)

Draw the structure with twist order parameter.

version 0.1.2

## Requirements

* genice<2.0
* genice-svg>=0.4
* twist-op>=0.2
* sklearn

## Installation from PyPI

    % pip install genice_twist

## Manual Installation

### System-wide installation

    % make install

### Private Installation

Copy the files in formats/ into your local formats/ folder.

## Usage

    
    usage:
        genice II -f twist[options:separated:by:colons] > file
    
    Output the twist values for all the hydrogen-bonded pairs.
    
    options:
        png      Draw the hydrogen bonds with a rainbow palette according to the twist value in PNG format.
        png:CM   Draw the hydrogen bonds with color-mixing scheme in PNG format.
        png:DB   Draw the hydrogen bonds with decision-boundary coloring scheme in PNG format.
        png:SB   Draw the hydrogen bonds with simple boundary coloring scheme in PNG format.
        svg      Draw the hydrogen bonds with a rainbow palette according to the twist value in SVG format.
        svg:CM   Draw the hydrogen bonds with color-mixing scheme in SVG format.
        svg:DB   Draw the hydrogen bonds with decision-boundary coloring scheme in SVG format.
        svg:SB   Draw the hydrogen bonds with simple boundary coloring scheme in SVG format.
        yaplot   Draw the hydrogen bonds with a rainbow palette according to the twist value in YaPlot format.
        shadow   Draw shadows to the atoms (PNG and SVG)
        Ih=filename.twhist   Specify the (two-dimensional) histogram of twist parameter in pure ice Ih.
        Ic=filename.twhist   Specify the (two-dimensional) histogram of twist parameter in pure ice Ic.
        LDL=filename.twhist  Specify the (two-dimensional) histogram of twist parameter in pure LDL.
        HDL=filename.twhist  Specify the (two-dimensional) histogram of twist parameter in pure HDL.
        rotatex=30   Rotate the picture (SVG and PNG)
        rotatey=30   Rotate the picture (SVG and PNG)
        rotatez=30   Rotate the picture (SVG and PNG)

## Auxiliary Files

* IhST2.twhist
* IcST2.twhist
* LDLST2.twhist
* HDLST2.twhist

They are two-dimensional (i.e. real and imaginary) histograms of twist order parameter of ST2 water at 235 K, 0.98 g cm$^{-3}$. You have to prepare the appropriate histogram if you want to apply CM (color-mixing) or DB (decision-boundary) coloring scheme to other water models or in other conditions.

## Test in place

    % make test

## Reference

* Matsumoto, M., Yagasaki, T. & Tanaka, H. A Bayesian approach for identification of ice Ih, ice Ic, high density, and low density liquid water with a torsional order parameter. J. Chem. Phys. 150, 214504 (2019).
