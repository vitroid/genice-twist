# coding: utf-8
"""
Draw the structure with twist order parameter.

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
"""


desc = { "ref": { "MYT2019": 'Matsumoto, M., Yagasaki, T. & Tanaka, H. A Bayesian approach for identification of ice Ih, ice Ic, high density, and low density liquid water with a torsional order parameter. J. Chem. Phys. 150, 214504 (2019).'},
         "brief": "Twist order parameter.",
         "usage": __doc__,
         }


from math import atan2, sin, cos, pi
import cmath
import colorsys
from logging import getLogger
import io
import sys

import numpy as np
import networkx as nx
import yaplotlib as yp

from genice2_svg.formats.svg import draw_cell
from genice2_svg.render_png import Render as pRender
from genice2_svg.render_svg import Render as sRender
import twist_op as top
import genice2.formats
from genice2.decorators import timeit, banner


class Twist():
    def __init__(self, graph, relcoord, cell):
        self.graph = graph
        self.relcoord = relcoord
        self.cell = cell


    def iter(self):
        logger = getLogger()
        vecs = np.zeros([len(self.graph.edges()),3])
        cent = dict()
        edges = self.graph.edges()
        for i, edge in enumerate(edges):
            a,b = edge
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            ad = np.dot( d, self.cell )
            ad /= np.linalg.norm(ad)
            vecs[i] = ad
            c = self.relcoord[a] + d/2
            ac = np.dot( c, self.cell )
            cent[edge] = ac
        for edge, tw in top.twist_iter(edges, vecs):
            yield edge, cent[edge], tw


    def serialize(self, tag):
        s = "# i, j, center, of, bond, bond-twist\n"
        s += "{0}\n".format(tag)
        s += self.relcoord.shape[0].__str__()+"\n"
        for edge,center,twist in self.iter():
            a,b = edge
            s += "{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f}\n".format(a,b,*center*10, twist.real, twist.imag)
        s += "-1 -1 0 0 0 0 0\n"
        return s

    def yaplot(self):
        logger = getLogger()
        maxcir = 16
        maxrad = 6
        # 16x6=96 color variations
        s = ""
        for cir in range(maxcir):
            for rad in range(maxrad):
                angle = cir*360/maxcir
                hue = ((angle - 60 + 360) / 360) % 1.0
                bri = rad / (maxrad-1)
                sat = sin(angle*pi/180)**2
                logger.debug((angle,sat,bri))
                r,g,b = colorsys.hsv_to_rgb(hue, sat, bri)
                n = cir*maxrad + rad + 3
                s += yp.SetPalette(n,int(r*255),int(g*255),int(b*255))
        for pair,center,twist in self.iter():
            if twist == 0:
                # not an appropriate pair
                continue
            a,b = pair
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], self.cell)
            bpos = apos + np.dot(d, self.cell)
            cosine = twist.real
            sine   = twist.imag
            angle = atan2(sine,cosine) * 180 / pi
            if angle < 0:
                angle += 360
            cir = int(angle*maxcir/360+0.5)
            # rad is squared.
            rad = int(abs(twist)**2 * maxrad)
            if cir > maxcir-1:
                cir -= maxcir
            if rad > maxrad-1:
                rad = maxrad-1
            palette = cir*maxrad + rad + 3
            # logging.info((abs(twist),rad,cir,palette))
            s += yp.Color(palette)
            s += yp.Line(apos,bpos)
            s += "# angle {0} rad {1} cir {2} rad {3}\n".format(angle, abs(twist), hue, rad)
        return s


    def svg(self, rotmat, render=sRender, shadow=None):
        Rsphere = 0.04  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for pair,center,twist in self.iter():
            a,b = pair
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            o = dp / np.linalg.norm(dp)
            o *= RR

            #Color setting
            cosine = twist.real
            sine   = twist.imag
            angle = atan2(sine,cosine) * 180 / pi
            rad = sine**2 + cosine**2
            hue = ((angle - 60 + 360) / 360) % 1.0
            bri = rad
            sat = sin(angle*pi/180)**2
            r,g,b = colorsys.hsv_to_rgb(hue, sat, bri)
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(r*255),int(g*255),int(b*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=shadow,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))



    def svg2(self, rotmat, phasefiles, render=sRender, shadow=None):
        #
        # Twist order parameter to distinguish phases.
        #
        logger = getLogger()

        pxF = dict()
        phases = list(phasefiles)
        logger.info(phases)
        X = []
        for phase in phases:
            pxF[phase] = np.loadtxt(open(phasefiles[phase], "r"))
            X.append(pxF[phase].reshape([1600,]))
        X = np.array(X).T

        chirs = []
        chiis = []
        for pair,center,twist in self.iter():
            chirs.append(twist.real)
            chiis.append(twist.imag)

        #pX: p(chi)
        pX = np.histogram2d(chirs, chiis,
                            bins=(40,40),
                            range=[[-1.,1.],[-1.,1.]],
                            normed=True)
        Y = pX[0].reshape([1600,])
        from sklearn import linear_model
        reg = linear_model.LinearRegression(fit_intercept=False)
        reg.fit(X,Y)

        # pF : p(F)

        pF = dict()
        for i, phase in enumerate(phases):
            pF[phase] = reg.coef_[i]

        pFx = dict()
        for phase in phases:
            pFx[phase] = pxF[phase]*pF[phase] / pX[0]
        pFx["ice"] = pFx["1c"] + pFx["1h"]

        Rsphere = 0.04  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for pair,center,twist in self.iter():
            a,b = pair
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            o = dp / np.linalg.norm(dp)
            o *= RR
            center = apos + dp/2

            #Color setting
            bin = int(twist.real*19.999+20), int(twist.imag*19.999+20)
            if pX[0][bin] == 0.0:
                continue
            green = pFx["HDL"][bin]
            blue  = pFx["LDL"][bin]
            red   = pFx["ice"][bin]
            logger.debug((red,green,blue,bin,pX[0][bin]))
            if green < 0:
                green = 0
            if green > 1:
                green = 1
            if red < 0:
                red = 0
            if red > 1:
                red = 1
            if blue < 0:
                blue = 0
            if blue > 1:
                blue = 1
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=shadow,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))


    def svg3(self, rotmat, render=sRender, shadow=None):
        #
        # Twist order parameter to distinguish phases.
        # Simple criteria.
        #
        logger = getLogger()
        anglerange = 20

        Rsphere = 0.04  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for pair,center,twist in self.iter():
            a, b = pair
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            o = dp / np.linalg.norm(dp)
            o *= RR
            center = apos + dp/2

            #Color setting
            cosine, sine  = twist.real, twist.imag
            angle = atan2(sine,cosine) * 180 / pi
            radius = (cosine**2 + sine**2)**0.5
            red = 0
            green = 0
            blue = 0
            if radius > 0.65:
                if -0.85 < twist.real < +0.85:
                    blue = 1
                else:
                    red = 1
            else:
                green = 1
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=shadow,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))



    def svg4(self, rotmat, phasefiles, render=sRender, shadow=None):
        #
        # Twist order parameter to distinguish phases.
        # Color only when likelihood > 0.5.
        #
        logger = getLogger()
        logger.info("Calculating OP...")
        pxF = dict()
        phases = list(phasefiles)
        logger.info(phases)
        X = []
        for phase in phases:
            pxF[phase] = np.loadtxt(open(phasefiles[phase], "r"))
            X.append(pxF[phase].reshape([1600,]))
        X = np.array(X).T

        chirs = []
        chiis = []
        for pair,center,twist in self.iter():
            chirs.append(twist.real)
            chiis.append(twist.imag)

        #pX: p(chi)
        pX = np.histogram2d(chirs, chiis,
                            bins=(40,40),
                            range=[[-1.,1.],[-1.,1.]],
                            normed=True)
        Y = pX[0].reshape([1600,])
        from sklearn import linear_model
        reg = linear_model.LinearRegression(fit_intercept=False)
        reg.fit(X,Y)

        # pF : p(F)

        pF = dict()
        for i, phase in enumerate(phases):
            pF[phase] = reg.coef_[i]

        pFx = dict()
        for phase in phases:
            pFx[phase] = pxF[phase]*pF[phase] / pX[0]
        pFx["ice"] = pFx["1c"] + pFx["1h"]
        logger.info("Done.")

        Rsphere = 0.04  # nm
        Rcyl    = 0.02  # nm
        RR      = (Rsphere**2 - Rcyl**2)**0.5
        prims = []
        proj = np.dot(self.cell, rotmat)
        xmin, xmax, ymin, ymax = draw_cell(prims, proj)
        for pair,center,twist in self.iter():
            a,b = pair
            d = self.relcoord[b] - self.relcoord[a]
            d -= np.floor(d + 0.5)
            apos = np.dot(self.relcoord[a], proj)
            dp = np.dot(d, proj)
            bpos = apos + dp
            o = dp / np.linalg.norm(dp)
            o *= RR
            center = apos + dp/2

            #Color setting
            bin = int(twist.real*19.999+20), int(twist.imag*19.999+20)
            green = pFx["HDL"][bin]
            blue  = pFx["LDL"][bin]
            red   = pFx["ice"][bin]

            if green > 0.5:
                green = 1.0
                blue = 0
                red = 0
            elif blue > 0.5:
                blue = 1
                red = 0.0
                green = 0
            elif red > 0.5:
                red = 1
                blue = 0
                green = 0
            else:
                green = 1
                blue = 1
                red = 1
            colorcode = "#{0:02x}{1:02x}{2:02x}".format(int(red*255),int(green*255),int(blue*255))
            prims.append([center, "L", apos+o, bpos-o,Rcyl, {"fill":colorcode}])
        for v in self.relcoord:
            prims.append([np.dot(v, proj),"C",Rsphere, {"fill":"#fff"}]) #circle
        return render(prims, Rsphere, shadow=shadow,
                   topleft=np.array((xmin,ymin)),
                   size=(xmax-xmin, ymax-ymin))

class Format(genice2.formats.Format):


    def __init__(self, **kwargs):
        logger = getLogger()
        logger.info("Hook0: ArgParser.")
        self.type = "file" # output the raw twist values in @BTWC format
        self.scheme = None # rainbow coloring
        self.shadow = None
        self.rotmat = np.array([[1., 0, 0], [0, 1, 0], [0, 0, 1]])
        self.phases = {"1h": "IhST2.twhist",
                             "1c": "IcST2.twhist",
                             "LDL": "LDLST2.twhist",
                             "HDL": "HDLST2.twhist"}
        for key, value in kwargs.items():
            logger.info("Option with arguments: {0} := {1}".format(key,value))
            if key == "rotmat":
                value = re.search(r"\[([-0-9,.]+)\]", value).group(1)
                self.rotmat = np.array([float(x) for x in value.split(",")]).reshape(3,3)
            elif key == "rotatex":
                value = float(value)*pi/180
                cosx = cos(value)
                sinx = sin(value)
                R = np.array([[1, 0, 0], [0, cosx, sinx], [0,-sinx, cosx]])
                self.rotmat = np.dot(self.rotmat, R)
            elif key == "rotatey":
                value = float(value)*pi/180
                cosx = cos(value)
                sinx = sin(value)
                R = np.array([[cosx, 0, -sinx], [0, 1, 0], [sinx, 0, cosx]])
                self.rotmat = np.dot(self.rotmat, R)
            elif key == "rotatez":
                value = float(value)*pi/180
                cosx = cos(value)
                sinx = sin(value)
                R = np.array([[cosx, sinx, 0], [-sinx, cosx, 0], [0, 0, 1]])
                self.rotmat = np.dot(self.rotmat, R)
            elif value is not True:
                #may be the file names for phases. just record them.
                self.phases[key] = value
                if value == "":
                    del self.phases[key]
            else:
                a = key
                logger.info("Flags: {0}".format(a))
                if a == "shadow":
                    self.shadow = "#4441"
                elif a in ("svg", "png", "yaplot"):
                    self.type = a
                elif a in ("CM", "DB", "SB"):
                    self.scheme = a
                else:
                    assert False, "Wrong options."
        logger.info("Hook0: end.")


    def hooks(self):
        return {2:self.Hook2}

    @timeit
    @banner
    def Hook2(self, lattice):
        "Complex bond twists."
        logger = getLogger()
        cell = lattice.repcell

        positions = lattice.reppositions
        graph = nx.Graph(lattice.graph) #undirected

        twist = Twist(graph, positions, cell.mat)
        s = ""
        if self.type == "yaplot":
            s = twist.yaplot()
        elif self.type == "svg":
            if self.scheme == "CM":
                s = twist.svg2(self.rotmat, self.phases)
            elif self.scheme == "SB":
                s = twist.svg3(self.rotmat)
            elif self.scheme == "DB":
                s = twist.svg4(self.rotmat, self.phases)
            else:
                # rainbow
                s = twist.svg(self.rotmat)
        elif self.type == "png":
            if self.scheme == "CM":
                s = twist.svg2(self.rotmat, self.phases, render=pRender, shadow=self.shadow)
            elif self.scheme == "SB":
                s = twist.svg3(self.rotmat, render=pRender, shadow=self.shadow)
            elif self.scheme == "DB":
                s = twist.svg4(self.rotmat, self.phases, render=pRender, shadow=self.shadow)
            else:
                s = twist.svg(self.rotmat, render=pRender, shadow=self.shadow)
        else: # "file"
            s = cell.serialize_BOX9()
            s += twist.serialize("@BTWC")
        self.output = s
