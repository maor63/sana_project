from __future__ import print_function

import timeit

start = timeit.default_timer()
top = 100000
for x in range(top):
    stop = timeit.default_timer()
    print('\r {}/{} time {}sec'.format(x, top, stop-start))