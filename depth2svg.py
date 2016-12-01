#!/usr/bin/env python3

import sys
import math
from collections import defaultdict
from argparse import ArgumentParser as AP


cli = AP(description='reads GR files and generates SVG')
cli.add_argument('-s', '--stagger', help='create staggering charts (default: off)',
                 type=int, default=1)
cli.add_argument('-w', '--width', help='width of bars (default: 5)', type=int, default=5)
cli.add_argument('-H', '--height', help='maximum height of bars (default: 100.0)',
                 type=float, default=100.0)
cli.add_argument('-d', '--distance', help='distance between bars (default: 2)',
                 type=int, default=2)
cli.add_argument('-t', '--ticks', help='number of grid lines (default: 2)',
                 type=int, default=2)
args = cli.parse_args()
width = args.width
height = args.height
space = args.distance
stagger = args.stagger
ticks = args.ticks

# TODO: inlcude sequence from fasta, features from gff and ticks
# TODO: inluce scale and horizontal bars
# make width fixed based on letter size, scaling should only be done by height
# multiple profiles makes no sense, remove that.


def nice_ticks(value, ticks=2):
    exponent = math.floor(math.log(value, 10))
    floor = math.floor(value / 10 ** exponent)
    if floor % ticks == 0 or floor / ticks % 0.5 == 0:
        return [i * floor * 10 ** exponent / ticks for i in range(1, ticks + 1)]
    if floor == 1:
        exponent -= 1
        floor = 11
    return nice_ticks((floor-1) * 10 ** exponent, ticks)


def write_svg(plus, minus):
    pos = 0
    max_bars = max([max(plus), max(minus)])
    max_v = max([max(plus.values()), -min(minus.values())])
    stack = int(max_v/stagger)+1
    pos = 0
    if stack <= 0:
        return
    print('''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%d" height="%d">
    <defs> <style type="text/css"><![CDATA[ %s %s ]]></style> </defs>'''
          % (max_bars*(width+space), height*2+space, default, css))
    if stagger == 1 and ticks > 0:
        for t in nice_ticks(max_v, ticks):
            h = height / max_v * t
            print('<line x1="0" y1="%.2f" x2="%d" y2="%.2f" class="tick" />' %
                  (height-h, max_bars*(width+space), height-h))
            print('<line x1="0" y1="%.2f" x2="%d" y2="%.2f" class="tick" />' %
                  (height+h, max_bars*(width+space), height+h))
    for i in range(0, max_bars+1):
        for p in range(1, int((plus[i]-1)/stack)+1):
            print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="plus" '
                  'fill-opacity="%.2f" id="p%d" />'
                  % (pos, 0, width, height, float(1/stagger), pos))
        h = plus[i] % stack / stack
        print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="plus" '
              'fill-opacity="%.2f" id="pp%d" />'
              % (pos, height-(height*h), width, height*h, float(1/stagger), pos))
        for m in range(1, int((-minus[i]-1)/stack)+1):
            print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="minus" '
                  'fill-opacity="%.2f" id="m%d" />'
                  % (pos, height+space, width, height, float(1/stagger), pos))
        h = -minus[i] % stack / stack
        print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="minus" '
              'fill-opacity="%.2f" id="pm%d" />'
              % (pos, height+space, width, height*h, float(1/stagger), pos))
        pos += width+space
    print('</svg>')

default = '''
.tick {
    stroke-width: 1;
    stroke: #aaa;
}
'''
css = ''

last_chro = 'seq'
plus = defaultdict(int)
minus = defaultdict(int)
plus[0] = 0
minus[0] = 0
for l in sys.stdin:
    c = l.split()
    if len(c) == 3:
        chro, pos, v = c[0], int(c[1]), int(c[2])
    elif len(c) == 2:
        chro, pos, v = 'seq', int(c[0]), int(c[1])
    else:
        sys.stderr.write('Input format not supported! Expected '
                         '"position<TAB>value" or "chromosome<TAB>position<TAB>value".\n')
        exit(1)

    if last_chro != 'seq' and chro != last_chro:
        write_svg(plus, minus)
        last_chro = chro
        plus = defaultdict(int)
        minus = defaultdict(int)
        plus[0] = 0
        minus[0] = 0

    if v > 0:
        plus[pos] += v
    else:
        minus[pos] += v

write_svg(plus, minus)
