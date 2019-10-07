#!/usr/bin/env python3
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


def qual_print(c, q):
    if ord(q)-33 > 34:
        sys.stdout.write('\033[34m%c' % c)
    elif ord(q)-33 > 32:
        sys.stdout.write('\033[36m%c' % c)
    elif ord(q)-33 > 30:
        sys.stdout.write('\033[32m%c' % c)
    elif ord(q)-33 > 28:
        sys.stdout.write('\033[33m%c' % c)
    elif ord(q)-33 > 26:
        sys.stdout.write('\033[31m%c' % c)
    else:
        sys.stdout.write('\033[35m%c' % c)


sys.stderr.write('\033[44m > 34 \033[46m > 32 \033[42m > 30 \033[43m > 28 \033[41m > 26 \033[45m <=26\033[0m\n')
for l in sys.stdin:
    if l.startswith('@'):
        sys.stdout.write(l.rstrip())
        s = sys.stdin.readline().rstrip()
        sys.stdin.readline()
        q = sys.stdin.readline().rstrip()
        sys.stdout.write('\t\033[1m')
        for i in range(0, len(s)):
            qual_print(s[i], q[i])
        sys.stdout.write('\t\033[0m\n')
