library(grid)
vp <- viewport(x=0.5,y=0.5,width=0.9, height=0.9)
pushViewport(vp)
grid.rect(gp=gpar(lty="dashed"))
grid.circle(x=0.6, y=0.4, r=0.3)

stickperson <- function() {
    grid.circle(x=.5, y=.8, r=.1, gp=gpar(fill="turquoise"))
    grid.lines(c(.5,.5), c(.7,.2)) # vertical line for body
    grid.lines(c(.5,.7), c(.6,.7)) # right arm
    grid.lines(c(.5,.3), c(.6,.7)) # left arm
    grid.lines(c(.5,.65), c(.2,0)) # right leg
    grid.lines(c(.5,.35), c(.2,0)) # left leg
}

stickperson()

vp1 <- viewport(x=0.5, y=0.75, width=0.6, height=0.3)
pushViewport(vp1)
grid.circle(gp=gpar(col="blue"))
# plot the outline of vp1:
grid.rect()

upViewport()
grid.circle(gp=gpar(col="purple"))

pushViewport(viewport())
grid.lines(c(.05, .95), c(.95, .05))
grid.lines(c(.05, .95), c(.05, .95))
for (i in 1:100) {
vp <- viewport(h=.9, w=.9)
pushViewport(vp)
grid.rect()
}

for (i in 1:30) {
    vp <- viewport(h=.9, w=.9)
    pushViewport(vp)
    # person 1:
    if(i == 5) {
pushViewport(viewport(x=.8))
     stickperson()
     upViewport()
     }
# person 2:
if(i == 20) {
pushViewport(viewport(x=.2))
     stickperson()
     upViewport()
     }
    # person 3:
    if(i == 30) stickperson()
}


escape_prop <-
c(0.24, 0.28, 0.28, 0.33, 0.33, 0.32, 0.3, 0.21, 0.3, 0.28, 0.17,
0.27, 0.21, 0.18, 0.22, 0.21, 0.19, 0.17, 0.17, 0.15, 0.25, 0.19,
0.19, 0.22, 0.21, 0.18, 0.24, 0.23, 0.27, 0.16, 0.17, 0.22, 0.17,
0.25, 0.19, 0.25, 0.12, 0.17, 0.22, 0.22)
nfires <-
c(953, 620, 584, 839, 1415, 1180, 656, 408, 872, 965, 853,
1492, 951, 772, 1541, 1114, 479, 860, 1166, 1208, 657, 1140,
1223, 1275, 489, 932, 1096, 1378, 1033, 889, 1046, 818, 1213,
782, 962, 1666, 2017, 1689, 1885, 1435)
nfirescode <- nfires/max(nfires)
index <- (1:40)/41

pushViewport(viewport(width=.9, height=.9))

pushViewport(viewport(y=.75, width=.9, height=.9))
for (i in 1:40) {
	vp <- viewport(x=index[i],y=escape_prop[i], height=.03, width=.03)
	pushViewport(vp)
	grid.circle(r=sqrt(nfirescode[i]/pi))
    upViewport()
}

grid.xaxis(at=c(0,index[c(10,20,30,40)]), label=seq(1960,2000,10))
grid.yaxis(at=seq(0,.5,.1))
grid.text("Proportion of Escaped Fires", y=.6)


burningtree <- function() {
    grid.rect(x=.5, y=.2, width=.2, height=.4, gp=gpar(fill="grey", col=NA))
    grid.circle(x=.5, y=.5, r=.3, gp=gpar(fill="orange", col=NA))
    pushViewport(viewport(clip="on"))
    pushViewport(viewport(x=.5, y=0, angle=45))
    grid.rect(x=.5, y=.5, width=.2, height=.2, gp=gpar(fill="grey", col=NA))
    upViewport(2)
}

pushViewport(viewport())
burningtree()

pushViewport(viewport(width=.9, height=.9))
pushViewport(viewport(y=.75, width=.9, height=.9))
for (i in 1:40) {
	vp <- viewport(x=index[i],y=escape_prop[i], height=nfirescode[i]/10,
    width=.03)
	pushViewport(vp)
    burningtree()  # this replaces the grid.circle of Figure 7
    upViewport()
}
grid.yaxis(at=seq(0,.5,.1))
grid.xaxis(at=c(0,index[c(10,20,30,40)]), label=seq(1960,2000,10))
grid.text("Proportion of Escaped Fires", y=.6)

gr <- rectGrob(width=0.1,height=0.1, name="gr")

gr1 <- editGrob(gr, vp=viewport(x=0.2, y=0.6), name="gr1")
gr2 <- editGrob(gr, vp=viewport(x=0.7, y=0.75), name="gr2")
gr3 <- editGrob(gr, vp=viewport(x=0.5, y=0.4), name="gr3")

grid.draw(gr1)
grid.draw(gr2)
grid.draw(gr3)

gr1 <- grid.edit("gr1", vp= viewport(x=.2,y=.6, angle=30))
gr2 <- grid.edit("gr2", vp= viewport(x=.7,y=.75, angle=63))
gr3 <- grid.edit("gr3", vp= viewport(x=.5,y=.4, angle=72))

for (i in 1:1000*5){
	grid.edit("gr1", vp= viewport(x=.2,y=.6, angle=i))
	grid.edit("gr2", vp= viewport(x=.7,y=.75, angle=i*2))
	grid.edit("gr3", vp= viewport(x=.5,y=.4, angle=i*3))
}

pushViewport(vp=viewport())
b2 <- sqrt(1/cos(36*pi/180)^2-1)/2
b3 <- sin(72*pi/180)/(2*(1+cos(72*pi/180))) - (1-sin(72*pi/180))/2
triangle2 <- polygonGrob(c(0,.5,1),c(b3,b2+b3,b3), name="triangle2", gp=gpar(fill="yellow", col=0))
grid.draw(triangle2)

for (i in 0:2){
	pushViewport(vp=viewport(angle=72*i))
	grid.draw(triangle2)
	upViewport()
}

grid.rect(gp=gpar(fill="blue")) # this gives a blue background
#draw star 1:
pushViewport(vp=viewport(x=.2, y=.2, w=.25, h=.25, angle=40))
for (i in 0:2){
	pushViewport(vp=viewport(angle=72*i))
	grid.draw(triangle2)
	upViewport(1)
}
upViewport(1)
#draw star 2:
pushViewport(vp=viewport(x=.8, y=.8, w=.3, h=.3, angle=90))
for (i in 0:2){
	pushViewport(vp=viewport(angle=72*i))
	grid.draw(triangle2)
	upViewport(1)
}
upViewport(1)
#draw star 3:
pushViewport(vp=viewport(x=.7, y=.3, w=.2, h=.2, angle=130))
for (i in 0:2){
	pushViewport(vp=viewport(angle=72*i))
	grid.draw(triangle2)
	upViewport(1)
}
upViewport(1)
#draw star 4:
pushViewport(vp=viewport(x=.3, y=.7, w=.15, h=.15, angle=210))
for (i in 0:2){
	pushViewport(vp=viewport(angle=72*i))
	grid.draw(triangle2)
	upViewport(1)
}
upViewport(1)

stickpersonGrob <- gList(circleGrob(x=.5, y=.8, r=.1, gp=gpar(fill="yellow")),
	linesGrob(c(.5,.5), c(.7,.2)), linesGrob(c(.5,.7), c(.6,.7)),
	linesGrob(c(.5,.3), c(.6,.7)), linesGrob(c(.5,.65), c(.2,0)),
	linesGrob(c(.5,.35), c(.2,0)))
grid.draw(stickpersonGrob)

stickpersonG <- gTree(children= stickpersonGrob, name="stickperson")
grid.draw(stickpersonG) # stick person on the left
grid.edit("stickperson", vp=viewport(angle=80)) # stick person on the right, rotated
