#!/usr/bin/env python3

import sys
from collections import defaultdict
from argparse import ArgumentParser as AP


cli = AP(description='reads GR files and generates SVG')
cli.add_argument('-s', '--stagger', help='create staggering charts (default: off)',
                 type=int, default=1)
cli.add_argument('-w', '--width', help='width of bars (default: 5)', type=int, default=5)
cli.add_argument('-m', '--max_height', help='maximum height of bars (default: 100.0)',
                 type=float, default=100.0)
cli.add_argument('-d', '--distance', help='distance between bars (default: 2)',
                 type=int, default=2)
args = cli.parse_args()
width = args.width
height = args.max_height
space = args.distance
stagger = args.stagger


def write_svg(plus, minus):
    pos = 0
    max_bars = max([max(plus), max(minus)])
    max_v = max([max(plus.values()), -min(minus.values())])
    stack = int(max_v/stagger)
    pos = 0
    if stack <= 0:
        return
    for i in range(0, max_bars+1):
        for p in range(1, int((plus[i]-1)/stack)+1):
            print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="plus" fill-opacity="%.2f" id="p%d" />'
                  % (pos, 0, width, height, float(1/stagger), pos))
        h = plus[i] % stack / stack
        print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="plus" fill-opacity="%.2f" id="pp%d" />'
              % (pos, height-(height*h), width, height*h, float(1/stagger), pos))
        for m in range(1, int((-minus[i]-1)/stack)+1):
            print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="minus" fill-opacity="%.2f" id="m%d" />'
                  % (pos, height+space, width, height, float(1/stagger), pos))
        h = -minus[i] % stack / stack
        print('<rect x="%d" y="%.2f" width="%d" height="%.2f" class="minus" fill-opacity="%.2f" id="pm%d" />'
              % (pos, height+space, width, height*h, float(1/stagger), pos))
        pos += width+space

default = ''
css = ''
print('''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="100%%" height="100%%">
      <defs> <style type="text/css"><![CDATA[ %s %s ]]></style> </defs>''' % (default, css))

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
print('</svg>')
