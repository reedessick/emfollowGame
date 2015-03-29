description = """ a module containing classes for emfollow game """

#=================================================

import pygame

import numpy as np
import healpy as hp

import sys
sys.path.append("/home/reed/LIGO/BAYESburst/skymap_statistics/")
import stats

#=================================================

### colors
black = np.array((  0,   0,   0))
white = np.array((255, 255, 255))
grey  = 0.5*white
red   = np.array((255,   0,   0))
green = np.array((  0, 255,   0))
blue  = np.array((  0,   0, 255))

### factors of pi
pi = np.pi
pi_2 = pi*0.5
twopi = 2*pi

#=================================================
class zoom( object ):

    def __init__(self, center, radius=100, crosshairs=False):
        self.x, self.y = center
        self.radius = radius
        self.diameter = 2*radius
        self.crosshairs = crosshairs

    def display(self, screen, color=black, fill=False):
        if fill:
            pygame.draw.circle(screen, white, (self.x, self.y), self.radius, 0) ### fill in 
        pygame.draw.circle(screen, color, (self.x, self.y), self.radius, 1) ### outline
        if self.crosshairs:
            pygame.draw.aaline(screen, color, (self.x+self.radius, self.y), (self.x-self.radius, self.y)) 
            pygame.draw.aaline(screen, color, (self.x, self.y+self.radius), (self.x, self.y-self.radius))

class mollweide( object ):
    dtheta = 30
    dphi = 45
    err = 1e-3
    maxdepth = 10

    def __init__(self, center, width=100, height=60, graticule=False, continents=False):
        self.x, self.y = center
        self.width = width
        self.hight = height

        self.halfwidth = width*0.5
        self.halfheight = height*0.5

        self.graticule = graticule
        self.continents = continents

    def isover(self, pos):
        x, y = self.localpos(pos)
        if abs(y) > 1: ### out of bounds vertically
            return False
        else:
            aa = np.arcsin(y)
            twoaa = 2*aa
            phi = pi + pi*x / np.cos(aa)

            return abs(phi-pi) <= pi
           
    def localpos(self, pos):
        return (pos[0]-self.x)/self.halfwidth, (pos[1]-self.y)/self.halfheight
 
    def display(self, screen, color=grey):
        ### draw boundaries 
        pts = np.concatenate( ( np.transpose(self.ang2pos( np.linspace(0, pi, 501), 0 )), np.transpose(self.ang2pos( np.linspace(pi, 0, 501), twopi )) ) )
        pygame.draw.polygon(screen, color, pts, 1 )

        if self.graticule:
            for phi in np.arange(45, 360, 45)*pi/180:
                pts = np.transpose(self.ang2pos( np.linspace(0, pi, 501), phi ))
                pygame.draw.aalines(screen, color, False, pts)
#                pygame.draw.lines(screen, color, False, pts)
            for theta in np.arange(30, 180, 30)*pi/180:
                pygame.draw.aaline(screen, color, self.ang2pos(theta, 0), self.ang2pos(theta, twopi))
#                pygame.draw.line(screen, color, self.ang2pos(theta, 0), self.ang2pos(theta, twopi))

        if self.continents:
            print "WARNING: continents not yet implemented..."
            pass  

    def theta2lat(self, theta):
        return theta - pi_2

    def lat2theta(self, latitude):
        return latitude + pi_2

    def ang2pos(self, theta, phi):
        if not isinstance(theta, np.ndarray):
            theta = np.array([theta])
        if not isinstance(phi, np.ndarray):
            phi = np.array([phi])

        latitude = self.theta2lat(theta)
        aa = self.auxang( latitude )

        return np.array((self.x + self.halfwidth*(phi/pi-1)*np.cos(aa), self.y + self.halfheight*np.sin(aa))).astype(int)

    def pos2ang(self, x, y):
        if not isinstance(x, np.ndarray):
            x = np.array([x])
        if not isinstance(y, np.ndarray):
            y = np.array([y])

        x, y = self.localpos((x,y))

        aa = np.arcsin(y)
        twoaa = 2*aa

        phi = pi + pi*x / np.cos(aa)
        latitude = np.arcsin( (twoaa + np.sin(twoaa))/pi)

        return self.lat2theta( latitude), phi

    def auxang(self, latitude):
         aa = np.ones_like( latitude )

         truth = np.abs(latitude)==pi_2
         aa[truth] = latitude[truth]
         if np.sum(truth) < len(aa):
             ntruth = (1-truth).astype(bool)
             aa[ntruth] = self.__auxang( latitude[ntruth], pi*np.sin(latitude[ntruth]), 0)
         return aa
 
    def __auxang(self, aa, pisinlat, depth):
        if depth > self.maxdepth:
            raise ValueError("recursive depth in mollweide.__auxang larger than %d"%self.maxdepth)           
        twoaa = 2*aa
        daa = -(twoaa + np.sin(twoaa) - pisinlat)/(2 + 2*np.cos(twoaa))
        if np.all(np.abs(daa) < self.err):
            return aa + daa
        else:
            return self.__auxang(aa+daa, pisinlat, depth+1)


class skymap( object ):

    def __init__(self, mollproj, post, nest=False):
        self.mollweide = mollproj

        self.post = post
        self.Mpost = np.max(post)
        self.mpost = np.min(post)

        self.npix = len(post)
        self.nside = hp.npix2nside(self.npix)
        self.nest = nest

        self.pixels = self.initpixels()

    def initpixels(self):
        return [pixel(self.mollweide, pix, self.nside, nest=self.nest, value=self.post[pix]) for pix in xrange(self.npix)]

    def display(self, screen, mollcolor=grey):
        for pix in self.pixels:
            pix.display(screen, self.colormap(pix.value)) ### FIX ME. get color from pix.value???
        
        self.mollweide.display( screen, color=mollcolor ) ### draw this on top of pixels...

    def colormap(self, value):
        return white*(self.Mpost - value)/(self.Mpost - self.mpost)

    def zoom_pixels(self, pos, radius=10*pi/180, nest=False):
        x, y = pos
        theta, phi = self.mollweide.pos2ang(x, y)
        good_pixels = self.tree_search(theta, phi, range(12), radius=radius, nside=1)
        if nest:
            return good_pixels
        else:
            return [(hp.nest2ring(self.nside, pixNo), angs) for pixNo, angs in good_pixels]

    def tree_search(self, theta, phi, pixNos, radius=10*pi/180, nside=1):
        if nside==self.nside:
            good_pix = []
            for pixNo in pixNos:
                t, p = hp.pix2ang(nside, pixNo, nest=True)
                d_theta = np.arccos(stats.cos_dtheta(t, p, theta, phi))
                if d_theta < radius:
#                    good_pix.append( ( pixNo, self.rotate(t, p, theta, phi) ) )
                    good_pix.append( ( pixNo, [self.rotate(_t, _p, theta, phi) for _t, _p in np.transpose(hp.vec2ang(hp.boundaries(nside, pixNo, nest=True, step=10).transpose())) ] ) )

            return good_pix
        else:
            good_pix = []
            angres = 2*hp.nside2pixarea(nside)**0.5
            for pixNo in pixNos:
                t, p = hp.pix2ang(nside, pixNo, nest=True)
                d_theta = np.arccos(stats.cos_dtheta(t, p, theta, phi))
                if d_theta < radius + angres:
                    pixNo4 = pixNo*4
                    good_pix += self.tree_search(theta, phi, [pixNo4, pixNo4+1, pixNo4+2, pixNo4+3], radius=radius, nside=nside*2)
            return good_pix

    def rotate(self, t, p, t0, p0):
        d_theta = np.arccos(stats.cos_dtheta(t, p, t0, p0))
        x = np.sin(t)*np.sin(p-p0)
        y = np.cos(t)*np.sin(t0) - np.sin(t)*np.cos(t0)*np.cos(p-p0)

        return d_theta, np.array((x,y))/(x**2+y**2)**0.5

class pixel( object ):
    center = None
    angs = None
    pts = None

    def __init__(self, mollproj, pix, nside, nest=False, value=0, dphi=0, selected=0):
        self.mollweide = mollproj

        self.value = value

        self.pix = pix
        self.nside = nside
        self.angres = hp.nside2pixarea(nside)**0.5
        self.nest = nest

        self.neighbors = [pixNo for pixNo in hp.get_all_neighbours(nside, pix, nest=nest) if pixNo+1 > 0]

        self.shift(dphi)

        self.selected = selected

    def isover(self, pos):
        if not self.mollweide.isover(pos):
            return False
        t, p = self.mollweide.pos2ang(pos[0], pos[1])

        ### hacky for now... check that this is closer than all it's neighbors
        this_dtheta = np.arccos(stats.cos_dtheta(self.center[0], self.center[1], t, p))
        neighbors = np.arccos(stats.cos_dtheta(self.neighbor_angs[:,0], self.neighbor_angs[:,1], t, p))
        return np.all(this_dtheta < neighbors)

        return np.arccos(stats.cos_dtheta(self.center[0], self.center[1], t, p)) < self.angres

    def shift(self, dphi):
        self.center = np.array(hp.pix2ang(self.nside, self.pix))
        self.center[1] += dphi

        self.neighbor_angs = np.transpose(hp.pix2ang(self.nside, self.neighbors))
        self.neighbor_angs[:,1] += dphi

        self.pts = self.getpts()

    def getangs(self, dphi=0):
        theta, phi = hp.vec2ang(hp.boundaries(self.nside, self.pix, nest=self.nest).transpose())
        #theta, phi = hp.vec2ang(hp.boundaries(self.nside, self.pix, nest=self.nest, step=10).transpose())
        return np.array((theta, phi+dphi))

    def getpts(self):
        if self.angs==None:
            self.angs = self.getangs()
        return np.array(self.mollweide.ang2pos(self.angs[0], self.angs[1])).transpose().astype(int)

    def display(self, screen, color=red):
        if self.pts==None:
            self.pts = self.getpts()

        pygame.draw.polygon(screen, color, self.pts)
        self.outline(screen, color=color)

    def outline(self, screen, color=red):
        if self.selected:
            pygame.draw.polygon(screen, color, self.pts, 1 )
        else:
            pygame.draw.polygon(screen, white, self.pts, 1 )
