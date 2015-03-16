#! /usr/bin/env python

import time, sys

if sys.platform == "win32":
    timer = time.clock
else:
    timer = time.time

class PTimer:
    t0 = 0
    t1 = 0

    def start(self):
        self.t0 = timer()

    def finish(self):
        self.t1 = timer()

    def seconds(self):
        return float(self.t1 - self.t0)






