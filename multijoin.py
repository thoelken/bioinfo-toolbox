#!/usr/bin/env python3
from argparse import ArgumentParser as AP
import sys
import locale

cli = AP('Multijoin', description='joins sorted files on one column')
cli.add_argument('-c', '--column', type=int, default=1, metavar='N',
                 help='column number to be joined on')
cli.add_argument('-d', '--delimiter', default='\t', metavar='STR',
                 help='delimiter for columns [default: "TAB" for \\t]')
cli.add_argument('-e', '--empty', default='', metavar='STR', help='fill empty fields')
cli.add_argument('files', nargs='+', help='files to join')
args = cli.parse_args()
args.column -= 1
if args.delimiter.lower() == 'tab':
    args.delimiter = '\t'

files = [open(f) for f in args.files]
lines = [f.readline().rstrip().split(args.delimiter) for f in files]
lines = [(l.pop(args.column), l) for l in lines]
keys, values = map(list, zip(*lines))
last = ''
while files:
    ordered = sorted(keys, key=locale.strxfrm)
    for c in ordered:
        current = c
        if c:
            break
    if not current:
        break
    if sorted([current, last], key=locale.strxfrm)[0] == current:
        sys.stderr.write('something is fishy with the order of %s and %s\n' % (current, last))
    last = current
    sys.stdout.write(current)
    for i in range(len(files)):
        sys.stdout.write(args.delimiter)
        if keys[i] == current:
            sys.stdout.write(args.delimiter.join(values[i]))
            l = files[i].readline()
            if not l:
                keys[i] = ''
            else:
                l = l.rstrip().split(args.delimiter)
                keys[i], values[i] = (l.pop(args.column), l)
        else:
            sys.stdout.write(args.delimiter.join([args.empty] * len(values[i])))
    print('')

for f in files:
    f.close()
