import sys, argparse
import math, os
#import numpy as np

class Helix(object):
    """
    Class to handle the three kinds of helices.
    Containing draw commands and graphical parameters.

    The simple bezier for the 1st-quarter sine in block size (1,1) is:
    0     0
    1/pi  0.5
    2/pi  1.0
    1.0   1.0

    Use SVG definitions to cheat with color name conversion.
    Available names are: ref
    """

    def __init__(self, nres, w, h):
        self.nres=nres
        self.minRes=3
        self.hMult=1.0
        self.hand=1.0
        self.rWidth=w
        self.rHeight=h

    def _qsine(self, xinit, xfact=1.0, yfact=1.0):
        """Slope first."""
        c0=(xfact*self.rWidth/math.pi, yfact*self.rHeight/2)
        c1=(xfact*2*self.rWidth/math.pi, yfact*self.rHeight)
        c2=(xfact*self.rWidth, yfact*self.rHeight)
        xW=xfact*self.rWidth
        if  xfact > 0.0:
            ostr='d="M %g,%g ' % (xinit[0],xinit[1])
        elif xfact < 0.0:
            ostr='d="M %g,%g ' % (xinit[0]+self.rWidth*(1.0+math.fabs(xfact)),xinit[1])
        ostr=ostr+'c %g,%g %g,%g %g,%g h %g c %g,%g %g,%g %g,%g z"' % (
                c0[0], c0[1], c1[0], c1[1], c2[0], c2[1], self.rWidth,
                c1[0]-xW, 0.0, c0[0]-xW, -1*c0[1], -1*c2[0], -1*c2[1] )
        return ostr

    def _hsine(self, xinit, xfact=1.0, yfact=1.0):
        ostr='d="M %g,%g ' % (xinit[0],xinit[1]-self.rHeight*yfact)
        c0=(xfact*self.rWidth/math.pi, yfact*self.rHeight/2.0)
        c1=(xfact*2.0*self.rWidth/math.pi, yfact*self.rHeight)
        c2=(xfact*self.rWidth, yfact*self.rHeight)
        xW=xfact*self.rWidth
        ostr=ostr+'c %g,%g %g,%g %g,%g %g,%g, %g,%g %g,%g h %g c %g,%g %g,%g %g,%g %g,%g, %g,%g %g,%g z"' % (
             xW-c1[0], 0.0, xW-c0[0], c0[1], c2[0], c2[1],
             c0[0], c0[1], c1[0], c1[1], c2[0], c2[1], self.rWidth,
             c1[0]-xW, 0.0, c0[0]-xW, -1.0*c0[1], -1.0*c2[0], -1.0*c2[1],
             -1.0*c0[0], -1.0*c0[1], -1.0*c1[0], -1.0*c1[1], -1.0*c2[0], -1.0*c2[1])
        return ostr

    def _crescent(self, xinit, xfact=1.0, yfact=1.0, yscale=0.6):
        ostr='d="M %g,%g ' % (xinit[0],xinit[1])
        c0=(xfact*self.rWidth/math.pi, yfact*self.rHeight/2.0)
        c1=(xfact*2.0*self.rWidth/math.pi, yfact*self.rHeight)
        c2=(xfact*self.rWidth, yfact*self.rHeight)
        xW=xfact*self.rWidth
        ostr=ostr+'c %g,%g %g,%g %g,%g %g,%g %g,%g %g,%g c %g,%g %g,%g %g,%g %g,%g %g,%g, %g,%g z"' % (
                c0[0], c0[1], c1[0], c1[1], c2[0], c2[1],
                xW-c1[0], 0.0, xW-c0[0], -1*c0[1], c2[0], -1*c2[1],
                -1*c0[0], yscale*c0[1], -1*c1[0], yscale*c1[1], -1*c2[0], yscale*c2[1],
                c1[0]-xW, 0.0, c0[0]-xW, -1*yscale*c0[1], -1*c2[0], -1*yscale*c2[1])
        return ostr

# = = = = = Temporary draw operations.
def print_define_linear_gradient(fp, name, box, stops):
    ntot=len(stops)
    print >> fp, '  <linearGradient id="%s" x1="%g" x2="%g" y1="%g" y2="%g">' % (name, box[0][0], box[1][0], box[0][1], box[1][1])
    for i in range(ntot):
        print >> fp, '    <stop offset="%g%%" stop-color="%s"/>' % (stops[i][0], stops[i][1])
    print >> fp, '  </linearGradient>'

def print_define_radial_gradient(fp, name, stops):
    ntot=len(stops)
    print >> fp, '  <radialGradient id="%s">' % (name)
    for i in range(ntot):
        print >> fp, '    <stop offset="%g%%" stop-color="%s"/>' % (stops[i][0], stops[i][1])
    print >> fp, '  </radialGradient>'

def print_define_helix_gradients(fp, suff, col):
    """
    col is a three-member list, going from dark to light / back to front.
    """
    print_define_linear_gradient(fp, 'GradHelixUF_'+suff, [[0,0],[0,1]], [ [0,col[1]],[100,col[2]] ] )
    print_define_linear_gradient(fp, 'GradHelixLF_'+suff, [[0,0],[0,1]], [ [0,col[2]],[100,col[1]] ] )
    print_define_linear_gradient(fp, 'GradHelixUB_'+suff, [[0,0],[0,1]], [ [0,col[1]],[100,col[0]] ] )
    print_define_linear_gradient(fp, 'GradHelixLB_'+suff, [[0,0],[0,1]], [ [0,col[0]],[100,col[1]] ] )
    print_define_linear_gradient(fp, 'GradHelixFF_'+suff, [[0,0],[0,1]], [[0,col[1]],[50,col[2]],[100,col[1]]])
    print_define_linear_gradient(fp, 'GradHelixFB_'+suff, [[0,0],[0,1]], [[0,col[1]],[45,col[0]],[55,col[0]],[100,col[1]]])

def print_defines(fp):
    print >> fp, '  <defs>'
    print_define_linear_gradient(fp, 'GradCoil',  [[0,0],[0,1]], [[0,'grey'],[50,'white'],[100,'grey']])
    print_define_radial_gradient(fp, 'GradTurn', [[0,'white'],[40,'lightcyan'],[100,'cyan']])
    print_define_linear_gradient(fp, 'GradSheet', [[0,0],[0,1]], [[0,'gold'],[50,'moccasin'],[100,'gold']])
    print_define_linear_gradient(fp, 'GradBridge', [[0,0],[0,1]], [[0,'darkkhaki'],[50,'tan'],[100,'darkkhaki']])
    print_define_helix_gradients(fp, 'H', ('darkred', 'red', 'orange') )
    print_define_helix_gradients(fp, 'G', ('darkmagenta', 'plum', 'pink') )
    print_define_helix_gradients(fp, 'I', ('indigo', 'purple', 'violet') )
    print >> fp, '  </defs>'

def print_svgheader(fp, width, height):
    print >> fp, '<?xml version="1.0" standalone="no"?>'
    print >> fp, '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"'
    print >> fp, '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'
    print >> fp, '<svg width="%g" height="%g" version="1.1" xmlns="http://www.w3.org/2000/svg">' % (width, height)

def print_svgfooter(fp):
    print >> fp, '</svg>'

def draw_coil(fp, xinit, length, thickness=15, linewidth=5):
    w=length ; th=thickness
    box=[ [ xinit[0],xinit[1]-th],[xinit[0]+w, xinit[1]+th] ]
    print >> fp, '  <path d="M %g,%g h %g v %g h %g z"' % ( xinit[0], xinit[1]-th, w, 2*th, -w)
    print >> fp, '    style="fill:url(#GradCoil)" />'
    print >> fp, '  <path d="M %g,%g h %g"' % ( xinit[0], xinit[1]-th, w )
    print >> fp, '    style="stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    print >> fp, '  <path d="M %g,%g h %g"' % ( xinit[0], xinit[1]+th, w )
    print >> fp, '    style="stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw

def draw_sheet(fp, x_init, width=100, nres=4, thickness=15, linewidth=5):
    w=width ;  th=thickness ; gap=0.3 ; omgap=1.0-gap
    draw_coil(fp, x_init, gap*w, th)
    draw_coil(fp, [x_init[0]+w*(nres-gap), x_init[1]], gap*w, th)
    print >> fp, '  <path d="M %g,%g l %g,%g h %g v %g l %g,%g ' % ( x_init[0]+gap*w, x_init[1]-th, omgap*w, -2.0*th, w*(nres-2), -2.0*th, omgap*w, 4.0*th)
    print >> fp, '    v %g l %g,%g v %g h %g l %g,%g z"' % ( 2.0*th, -omgap*w, 4.0*th, -2.0*th, w*(2-nres), -1*omgap*w, -2.0*th)
    print >> fp, '    style="fill:url(#GradSheet);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw

def draw_bridge(fp, x_init, width=100, thickness=15, linewidth=5):
    w=width ; th=thickness ; gap=0.2 ; omgap=1.0-gap ; hmgap=0.5-gap
    draw_coil(fp, x_init, gap*w, th)
    print >> fp, '  <path d="M %g,%g l %g,%g v %g l %g,%g ' % ( x_init[0]+gap*w, x_init[1]-th, hmgap*w, -2.0*th, -2.0*th, hmgap*w, 4.0*th)
    print >> fp, '    v %g l %g,%g v %g l %g,%g z"' % ( 2.0*th, -1*hmgap*w, 4.0*th, -2.0*th, -1*hmgap*w, -2.0*th)
    print >> fp, '    style="fill:url(#GradBridge);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    draw_coil(fp, [x_init[0]+w*omgap, x_init[1]], gap*w, th)

def draw_turn(fp, x_init, width=100, nres=4, thickness=15, linewidth=5):
    w=width
    th=thickness
    gap=0.2 ; dgap=2.0*gap
    print >> fp, '  <path d="M %g,%g h %g v %g h %g"' % (x_init[0], x_init[1]-th, gap*w, 2*th, -1*gap*w)
    print >> fp, '    style="fill:url(#GradCoil);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    print >> fp, '  <path d="M %g,%g a %g,%g %i %i %i %g,%g v %g a %g,%g %i %i %i %g,%g z"' % (
                  x_init[0]+gap*w, x_init[1]-th, 0.75*w*(nres-dgap), 2.5*th, 0, 0, 1, w*(nres-gap*2), 0.0,
                  2*th, 0.75*w*(nres-dgap), 2.5*th, 0, 0, 1, w*(dgap-nres), 0.0)
    print >> fp, 'style="fill:url(#GradTurn);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    print >> fp, '  <path d="M %g,%g h %g v %g h %g"' % (x_init[0]+w*nres, x_init[1]-th, -1*gap*w, 2*th, gap*w)
    print >> fp, '    style="fill:url(#GradCoil);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw

# Define only right handed helix for now.
def draw_helix(fp, x_init, width=100, nres=4, thickness=15, htype='H', linewidth=5):
    """
    This is the most complicated draw-call in the script. First, the helix must be analysed so as to scale it to the nearest half-turn.
    Then we print in the following order:
    1) Connecting rods to the surrounding coils.
    2) Back-halves of helices.
    3) Front-halves of helices,
    4) Thickness crescents to give it a 3D-look.
    An X-scale is needed to squeexe or stretch the helix without too much effort.
    """
    w=width ; th=thickness
    yscale=0.6
    if htype=='H':
        res_per_turn=3.6 ; hRad=5.0*th
    elif htype=='I':
        res_per_turn=4.1 ; hRad=5.5*th
    elif htype=='G':
        res_per_turn=3.1 ; hRad=4.5*th

    halfturns=int(round(2.0*(nres-1)/res_per_turn))

    xfact=0.5*(nres-1)/halfturns
    motif=Helix(nres, w, hRad )
    #print "Debug number of half-turns, xfact, w:", halfturns, xfact, w

    #Layer-0, Connecting tubes
    print >> fp, '  <path d="M %g,%g h %g l %g,%g h %g"' % (x_init[0], x_init[1]-th, w-th, 2*th, 2*th, -1*(w+th))
    print >> fp, '    style="fill:url(#GradCoil);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    if halfturns%2==0:
        print >> fp, '  <path d="M %g,%g h %g l %g,%g h %g"' % (x_init[0]+nres*w, x_init[1]-th, -1*(th+w), 2*th, 2*th, w-th)
        print >> fp, '    style="fill:url(#GradCoil);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw
    else:
        print >> fp, '  <path d="M %g,%g h %g l %g,%g h %g"' % (x_init[0]+nres*w, x_init[1]-th, th-w, -2*th, 2*th, w+th)
        print >> fp, '    style="fill:url(#GradCoil);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % lw

    #Layer-1, full back-twist
    for i in range(1, halfturns, 2):
        yfact=-1.0*yscale
        x0=x_init[0]+(2*i-1)*xfact*w
        #print "...Debug back:", i, x0
        ostr=motif._hsine( [ x0, x_init[1]], xfact=xfact, yfact=yfact)
        print >> fp, '<path %s' % ostr
        print >> fp, 'style="fill:url(#GradHelixFB_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( htype, lw )

    # Last twist for half turns.
    if halfturns%2==1:
        yfact=1.0*yscale
        x0=x_init[0]+(2*halfturns-2)*xfact*w
        ostr=motif._qsine( [ x0, x_init[1]], xfact=-1*xfact, yfact=yfact)
        print >> fp, '<path %s' % ostr
        print >> fp, 'style="fill:url(#GradHelixLB_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( htype , lw )

    # Layer=2, full front-twist
    for i in range(3, halfturns+1, 2):
        yfact=1.0
        x0=x_init[0]+(2*i-3)*xfact*w
        # print "...Debug front:", i, x0
        ostr=motif._hsine( [ x0, x_init[1]], xfact=xfact, yfact=yfact)
        print >> fp, '<path %s' % ostr
        print >> fp, 'style="fill:url(#GradHelixFF_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( htype, lw )

    # First twist.
    ostr=motif._qsine( x_init, xfact=xfact)
    print >> fp, '<path %s' % ostr
    print >> fp, 'style="fill:url(#GradHelixLF_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( htype, lw )
    # Last twist for full turns.
    if halfturns%2==0:
        x0=x_init[0]+w*(nres-2-xfact)
        ostr=motif._qsine( [ x0, x_init[1]], xfact=-1*xfact, yfact=-1.0)
        print >> fp, '<path %s' % ostr
        print >> fp, 'style="fill:url(#GradHelixUF_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( htype, lw )

    # Layer-3 thickness mods.
    for i in range(halfturns):
        if i%2==1:
            yfact=-1.0 ; lett='U'
            x0=x_init[0]+w*2*i*xfact
        else:
            x0=x_init[0]+w*(1+2*i*xfact)
            yfact=1.0 ; lett='L'
        ostr=motif._crescent( [ x0, x_init[1]], xfact=xfact, yfact=yfact, yscale=yscale)
        print >> fp, '<path %s' % ostr
        print >> fp, 'style="fill:url(#GradHelix%sF_%s);stroke:#000000;stroke-width:%g;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />' % ( lett, htype, lw )

# = = = = = Mundane operations like file I/O.

def read_dssp(fn):
    bData=False
    slist=''
    for l in open(fn):
        if len(l)<10:
            continue
        if not bData:
            if "#" in l:
                bData=True
            continue
        else:
            c=l[11] ; rn=l[13] ; ss=l[16]
            if rn != '!':
                ri=float(l[5:10])
            else:
                ss='!'
            if ss==" ":
                ss='C'
            slist+=ss
    return slist

def read_fastalike(fn):
    slist=''
    for l in open(fn):
        if l[0]=='>' or l=='':
            continue
        else:
            slist+=l[:-1]
    return slist

def interpret_ss_string(ss):
    ntot=len(ss)
    cprev=''
    elements=[]
    lengths=[]
    l=1 ; elements.append(ss[0])
    for i in range(1,ntot):
        if cprev!=ss[i]:
            lengths.append(l) ; l=1
            elements.append(ss[i])
        else:
            l+=1
        cprev=ss[i]
    lengths.append(l)
    return [ elements, lengths ]

# = = = = = = = = = = = = = = = =
# = = = = =  Main Program = = = =
# = = = = = = = = = = = = = = = =

scriptname=os.path.basename(__file__)
parser = argparse.ArgumentParser(description='Secondary-structure SVG Drawer. Takes as input one of several kinds of '
                                             'file formats, and extracts just the secondary-structure signature from this. '
                                             'Uses the DSSP-naming conventions as a basis for calculating, except "S" '
                                             'bends, which are resolved to coil.'
                                             'Also, if one desires PDF output, one way on linux is to run the following: '
                                             '"rsvg-convert -f pdf -o $i.pdf $i.svg" '
                                             'Copyright-2017 Poker Chen, until I release it on Github or something.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--dssp', type=str, dest='in_dssp', default='',
                    help='Input a DSSP file. This will search the file for a "#" character '
                         'indicating the start of the sequence data. Then reads column 17 which is its SS-assignment.')
parser.add_argument('-f', '--fasta', type=str, dest='in_fasta', default='',
                    help='Input a FASTA-like file, which the script will just scan and add characters from any line that does '
                         'not begin with ">".')
parser.add_argument('-o', '--out', type=str, dest='ofile', default='output.svg',
                    help='Output .svg file containing the scalable vector graphics.')
parser.add_argument('-s','--scale',type=float,default=1.0,
                    help='Overall scale of the output image, multiplying all indiviudal values below.')
parser.add_argument('-w','--width', type=float, default=100,
                    help='Unitless width of the structure taken up by each residue. Total width of the image will be Nres*w.')
parser.add_argument('-l','--height', type=float, default=200,
                    help='Unitless height of the output image.')
parser.add_argument('-t','--thickness', type=float, default=15,
                    help='Unitless thickness of the backbone-strand. This is used to determine the perpendicular size of'
                         'structure features. E.g., a helix will be ~5*th wide.' )
parser.add_argument('-lw', '--linewidth', type=float, default=5,
                    help='Width of borders around all graphical elements.')


args = parser.parse_args()
width=args.width*args.scale
height=args.height*args.scale
th=args.thickness*args.scale
lw=args.linewidth*args.scale

outfile=args.ofile
if args.in_dssp != '':
    ss = read_dssp(args.in_dssp)
elif args.in_fasta != '':
    ss = read_fastalike(args.in_fasta)
else:
    print '= = ERROR: This script needs an input file!'
    sys.exit(1)

print "= = Read input SS-sequence:"
print ss

structures = interpret_ss_string(ss)
ntot=len(ss)
nstruct=len(structures[0])

fp= open(outfile, 'w')
print_svgheader(fp, width*ntot, height)
print_defines(fp)
x_init=[0,0.5*height]

for i in range(nstruct):
    stype=structures[0][i]
    slen=structures[1][i]
    print '= = Debug:', x_init, stype, slen
    if stype=='C':
        draw_coil(fp, x_init, length=width*slen, thickness=th, linewidth=lw )
    elif stype=='G' or stype=='H' or stype=='I':
        draw_helix(fp, x_init, width=width, nres=slen, htype=stype, thickness=th, linewidth=lw)
    elif stype=='E':
        draw_sheet(fp, x_init, width, nres=slen, thickness=th, linewidth=lw)
    elif stype=='T':
        draw_turn(fp, x_init, width, nres=slen, thickness=th, linewidth=lw)
    elif stype=='B':
        draw_bridge(fp, x_init, width=width, thickness=th, linewidth=lw)
    elif stype=='!':
        #Assume DSSP-style chain-break
        print "= = NB: Chain break detected. Will skip."
    else:
        draw_coil(fp, x_init, width*slen, thickness=th, linewidth=lw)
    x_init[0]+=width*slen

print_svgfooter(fp)
fp.close()

print '= = = svg construction complete.'
