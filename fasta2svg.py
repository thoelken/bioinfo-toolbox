#!/usr/bin/env python3
import sys
from argparse import ArgumentParser as AP


cli = AP(description='''Create vector based alignment visualisation''')
cli.add_argument('-s', '--css', help='style sheet for colors and format')
cli.add_argument('-d', '--dump_default', help='print default CSS to screen', action='store_true')

default = '''
text {
    font-size:12px;
    font-style:normal;
    font-variant:normal;
    font-weight:normal;
    font-stretch:normal;
    text-align:start;
    letter-spacing:0px;
    word-spacing:0px;
    writing-mode:lr-tb;
    text-anchor:start;
    fill:#000000;
    fill-opacity:1;
    stroke:none;
    font-family:DejaVu Sans Mono;
}
text.name {
    font-size:10px;
    font-style:italic;
    text-align:end;
    font-family:DejaVu Sans;
    text-anchor:end;
}
rect {
    fill:#555;
    fill-opacity:0.5;
    fill-rule:nonzero;
    stroke:none;
}
rect.a {fill:#f77;}
rect.c {fill:#7f7;}
rect.g {fill:#77f;}
rect.t {fill:#ff3;}
'''
css = ''

args = cli.parse_args()
if args.dump_default:
    print(default)
    exit(1)
if args.css:
    with open(args.css, 'r') as f:
        css = f.read()

pos = 0
line = -1
names = []
boxes = []
letters = []
space = 5
for l in sys.stdin:
    if l.startswith('>'):
        line += 1
        names.append('''<text x="-10" y="%d" id="name%d" class="name">%s</text>'''
                     % (line*(space+12)+10, line, l.rstrip()[1:]))
        pos = 0
        continue
    if l == '' or l.startswith('#'):
        continue
    for c in l.rstrip():
        boxes.append('''<rect width="7" height="12" x="%d" y="%d" id="rect%d_%d" class="%s" />'''
                     % (pos*8, line*(space+12), pos, line, c.lower()))
        letters.append('''<text x="%d" y="%d" id="text%d_%d" class="%s">%s</text>'''
                       % (pos*8, line*(space+12)+10, pos, line, c.lower(), c))
        pos += 1

print('''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="100%%" height="100%%">
      <defs> <style type="text/css"><![CDATA[ %s %s ]]></style> </defs>''' % (default, css))

for n in names:
    print(n)
for b in boxes:
    print(b)
for l in letters:
    print(l)

print('</svg>')
