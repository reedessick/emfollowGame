description = """ play a game """

import pygame
pygame.init()

import utils

import sys
import time

#=================================================
### set up screen

screen_width = 1024
screen_height = 768
screen_size = (screen_width, screen_height)

screen = pygame.display.set_mode(screen_size)
screen.fill(utils.white)

clock = pygame.time.Clock()

#=================================================
### locations/sizes of elements

buff = 10
twobuff = 2*buff

o_width = 500
o_height = 300
o_center = buff + o_width/2, buff +o_height/2

#oz_radius = max(buff, 50)
#oz_center = buff + o_width/2, buff + o_height + buff + oz_radius

p_width = 500
p_height = 300
p_center = buff+ p_width/2 , buff + o_height + buff + p_height/2 

pz_radius = max(buff, 100)
pz_center = buff+p_width + buff + pz_radius, buff + o_height + buff + p_height/2

#========================
# text!
font = pygame.font.Font(None, 20)

### loading messages
#obs_loading_message = font.render("loading observatories...", 1, utils.black)
#screen.blit(obs_loading_message, o_center) 

map_loading_message = font.render("loading skymap", 1, utils.red)
screen.blit(map_loading_message, (p_center[0]-50, p_center[1])) 

obs_label = font.render('Observatories', 1, utils.black)
screen.blit(obs_label, (o_center[0]-o_width/2, o_center[1]-o_height/2))

map_label = font.render('Skymap', 1, utils.black)
screen.blit(map_label, (p_center[0]-p_width/2, p_center[1]-p_height/2))

zoom_label = font.render('Zoom', 1, utils.black)
screen.blit(zoom_label, (pz_center[0]-pz_radius, pz_center[1]-pz_radius))

pygame.display.flip()

#========================
### map of observatories
#o_zoom = utils.zoom(oz_center, radius=oz_radius, crosshairs=True)
o_moll = utils.mollweide( o_center, width=o_width, height=o_height, graticule=True, continents=True)

o_moll.display(screen)
pygame.display.flip()

#========================
### projection for skymaps
p_zoom = utils.zoom(pz_center, radius=pz_radius, crosshairs=True)
p_moll = utils.mollweide( p_center, width=p_width, height=p_height, graticule=True, continents=False)

fits = "sample16.fits"
post, header = utils.hp.read_map(fits, h=1)
header = dict(header)
p_map = utils.skymap( p_moll, post, nest=header['ORDERING']=='NEST' )

angres = utils.hp.nside2pixarea(utils.hp.npix2nside(len(post)))**0.5
p_zoomthr = 5*angres
zbuff = 0.05*pz_radius


p_map.display(screen)
#p_moll.display(screen)
pygame.display.flip()

#=================================================
### run the game
runme = True
framerate = 120

timer = 0
moved = False
stopthr = 10
mousevent=None
zoom_pix = []

while runme:

    clicked = False

    for event in pygame.event.get():

        if event.type == pygame.QUIT:
            runme = False

        elif hasattr(event, 'key'):
            if event.key == pygame.K_ESCAPE:
                runme = False

        elif event.type == pygame.MOUSEMOTION:
            moved = True
            timer = 0
            mousevent=event

        elif event.type == pygame.MOUSEBUTTONUP:
            clicked = True
            for pix in p_map.pixels:
                if pix.isover( event.pos ):
                    pix.selected = 1-pix.selected

    timer += clock.tick(framerate)

    ### update only if mouse has stopped...
    if moved and (timer >= stopthr):
        moved = False

        ### observatories!
        if o_moll.isover(mousevent.pos):
             pass
             ### need to detect if we've selected an observatory!
             ### handle all that jazz...

        ### posterior!
        pygame.draw.circle(screen, utils.white, (p_zoom.x, p_zoom.y), p_zoom.radius)
        if p_moll.isover(mousevent.pos):
            selected = []
            zoom_pix = []
            good_pixels = p_map.zoom_pixels(mousevent.pos, radius=p_zoomthr)
            for pixNo, angs in good_pixels:
                pts = []
                for x, y in [(x*dtheta*p_zoom.radius/p_zoomthr, -y*dtheta*p_zoom.radius/p_zoomthr) for dtheta, (x,y) in angs]:
                    r = (x**2 + y**2)**0.5
                    if r > p_zoom.radius:
                        x *= p_zoom.radius/r
                        y *= p_zoom.radius/r
                    pts.append( (p_zoom.x + x, p_zoom.y + y) ) 

                zoom_pix.append( (pixNo, pts) )

                pygame.draw.polygon(screen, p_map.colormap(p_map.pixels[pixNo].value), pts, 0)
                if p_map.pixels[pixNo].selected:
                    selected.append( pts )
                else:
                    pygame.draw.polygon(screen, utils.grey, pts, 1)
            for pts in selected:
                pygame.draw.polygon(screen, utils.red, pts, 1)

        p_zoom.display(screen, color=utils.grey, fill=False)

    if clicked:
       selected = []
       for pixNo, pts in zoom_pix:
           if p_map.pixels[pixNo].selected:
               selected.append( pts )
           else:
               pygame.draw.polygon(screen, utils.grey, pts, 1)

       for pts in selected:
           pygame.draw.polygon(screen, utils.red, pts, 1)
       
       for pix in p_map.pixels:
           pix.outline(screen, color=utils.red)

       p_moll.display(screen, color=utils.grey) 

    ### display!
    pygame.display.flip()

### quit the game!
pygame.quit()
#sys.exit(0)
