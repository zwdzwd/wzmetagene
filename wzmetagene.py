#!/usr/bin/env python
"""
Generate intervals for meta-gene plot.
Usage:
python ~/wzlib/pyutils/wzmetagene.py input_bed | les
"""

import argparse
import numpy as np

class Record(object):

    def __init__(self, fields, args):

        self.fields = fields
        self.chrm = fields[0]
        self.beg  = int(fields[1])
        self.end  = int(fields[2])
        
        if args.collapse:
            self.mid = (self.beg + self.end + 1) / 2
            self.beg = self.mid
            self.end = self.mid

        self.step1 = args.flankstep
        self.step2 = args.flankstep

    def sample_forward(self, args, index_func):

        if self.step2 < 0:
            return

        _window_beg = self.end
        for i in range(args.flanknumber):
            _window_end = _window_beg + self.step2

            if args.outer:
                window_end = int(_window_end)
                window_beg = window_end - 1
            elif args.middle:
                # trivial to awk, just for convenience
                window_mid = int((_window_beg + _window_end)/2.0)
                window_beg = window_mid-1
                window_end = window_mid
            else:
                window_beg = int(_window_beg)
                window_end = int(_window_end)

            index = index_func(i)
            if index < 0:
                if args.flankbygene:
                    reg = '(%d)-(%d)' % (-i-1, -i)
                elif args.flanktoneighbor:
                    if args.fold:
                        reg = '(%d)-(%d)%%' % (
                            float(-i-1)/args.flanknumber/2*100, float(-i)/args.flanknumber/2*100)
                    else:
                        reg = '(%d)-(%d)%%' % (
                            float(-i-1)/args.flanknumber*100, float(-i)/args.flanknumber*100)
                else:
                    reg = '(%d)-(%d)' % (self.end - window_end, self.end - window_beg)
            else:
                if args.flankbygene:
                    reg = '%d-%d' % (i, i+1)
                elif args.flanktoneighbor:
                    reg = '%d-%d%%' % (float(i)/args.flanknumber*100, float(i+1)/args.flanknumber*100)
                else:
                    reg = '%d-%d' % (window_beg - self.end, window_end - self.end)

            # in collapsed case, you can have index 0
            if index >= 0:
                area = 1
            else:
                area = -1

            if args.expansion > 0:
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion

            if window_beg > 0 and window_end > window_beg:
                print('%s\t%d\t%d\t%d\t%s\t%d\t%s' %
                      (self.chrm, window_beg, window_end,
                       index, reg, area, '\t'.join(self.fields)))
                
            _window_beg = _window_end
            

    def sample_backward(self, args, index_func):

        if self.step1 < 0:
            return

        _window_end = self.beg
        for i in range(args.flanknumber):
            _window_beg = _window_end - self.step1

            if args.outer:
                window_beg = int(_window_beg)
                window_end = window_beg + 1
            elif args.middle:
                # trivial to awk, just for convenience
                window_mid = int((_window_beg + _window_end)/2.0)
                window_beg = window_mid-1
                window_end = window_mid
            else:
                window_beg = int(_window_beg)
                window_end = int(_window_end)

            index = index_func(i)
            if index < 0:
                if args.flankbygene:
                    reg = '(%d)-(%d)' % (-i-1, -i)
                elif args.flanktoneighbor:
                    if args.fold:
                        reg = '(%d)-(%d)%%' % (
                            float(-i-1)/args.flanknumber/2*100, float(-i)/args.flanknumber/2*100)
                    else:
                        reg = '(%d)-(%d)%%' % (
                            float(-i-1)/args.flanknumber*100, float(-i)/args.flanknumber*100)
                else:
                    reg = '(%d)-(%d)' % (window_beg - self.beg, window_end - self.beg)
            else:
                if args.flankbygene:
                    reg = '%d-%d' % (i, i+1)
                elif args.flanktoneighbor:
                    if args.fold:
                        reg = '%d-%d%%' % (
                            float(i)/args.flanknumber/2*100, float(i+1)/args.flanknumber/2*100)
                    else:
                        reg = '%d-%d%%' % (
                            float(i)/args.flanknumber*100, float(i+1)/args.flanknumber*100)
                else:
                    reg = '%d-%d' % (self.beg - window_end, self.beg - window_beg)

            # in collapsed case, you can have index 0
            if index >= 0:
                area = 1
            else:
                area = -1

            if args.expansion > 0:
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion

            if window_beg > 0 and window_end > window_beg:
                print('%s\t%d\t%d\t%d\t%s\t%d\t%s' %
                      (self.chrm, window_beg, window_end,
                       index, reg, area, '\t'.join(self.fields)))

            _window_end = _window_beg

    def sample_internal(self, args, index_func):

        if args.outer:
            sentinels = list(np.linspace(self.beg+1, self.end, args.numinternal))
            for i in range(len(sentinels)):
                window_end = int(sentinels[i])
                window_beg = window_end - 1

                if args.expansion > 0:
                    window_beg = max(window_beg - args.expansion,0)
                    window_end = window_end + args.expansion

                index = index_func(i)
                if window_beg > 0 and window_end > window_beg:
                    print('%s\t%d\t%d\t%d\t%d-%d%%\t0\t%s' %
                            (self.chrm, window_beg, window_end, index,
                             float(index)/args.numinternal*100,
                             float(index+1)/args.numinternal*100, '\t'.join(self.fields)))
            return

        sentinels = list(np.linspace(self.beg, self.end, args.numinternal+1))
        for i in range(len(sentinels)-1):
            window_beg = int(sentinels[i])
            window_end = int(sentinels[i+1])
            if args.middle:
                window_mid = int((sentinels[i] + sentinels[i+1])/2)
                window_beg = window_mid-1
                window_end = window_mid

            if args.expansion > 0:
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion
                
            index = index_func(i)
            # 0 for internal
            if window_beg > 0 and window_end > window_beg:
                print('%s\t%d\t%d\t%d\t%d-%d%%\t0\t%s' %
                      (self.chrm, window_beg, window_end, index,
                       float(index)/args.numinternal*100,
                       float(index+1)/args.numinternal*100, '\t'.join(self.fields)))


def process_record(r, r0, r2):

    if args.flankbygene:
        r.step1 = float(r.end - r.beg + 1) / (args.numinternal)
        r.step2 = r.step1
    elif args.flanktoneighbor:
        if r0 is None or r.beg < r0.end:
            r.step1 = -1
        else:
            r.step1 = float(r.beg - r0.end + 1) / (args.flanknumber)
            if args.fold:
                r.step1 /= 2

        if r2 is None or r.end > r2.beg:
            r.step2 = -1
        else:
            r.step2 = float(r2.beg - r.end + 1) / (args.flanknumber)
            if args.fold:
                r.step2 /= 2

    if args.strand is None:
        strand = '+'
    else:
        strand = r.fields[args.strand-1]

    if strand == '+':
        if args.fold:
            r.sample_backward(args, lambda i: -i-1)
            r.sample_internal(args, lambda i: min(i, args.numinternal-1-i))
            r.sample_forward(args, lambda i: -i-1)
        else:
            r.sample_backward(args, lambda i: -i-1)
            r.sample_internal(args, lambda i: i)
            r.sample_forward(args, lambda i: args.numinternal+i)
    else:
        r.sample_forward(args, lambda i: -i-1)
        r.sample_internal(args, lambda i: args.numinternal-1-i)
        r.sample_backward(args, lambda i: args.numinternal+i)
        
    return

def main(args):

    if args.collapse:
        if args.outer: # inference: if only outer point is printed, then keep the collapsed middle
            args.numinternal = 1
        else:
            args.numinternal = 0

    r0 = None              # the previous record
    r  = None              # the current record
    r2 = None              # the next record
    for line in args.table:
        r2 = Record(line.strip('\n').split('\t'), args)

        if r is not None: # skip the first line
            if r.chrm == r2.chrm:
                process_record(r, r0, r2)
            else:               # chromosome switch
                process_record(r, r0, None)

        r0 = r
        r  = r2


    # the last record
    r2 = None
    process_record(r, r0, r2)

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate meta gene')
    parser.add_argument('table', help="Input bed file", type = argparse.FileType('r'), default='-')
    # report middle point for each sampled interval
    parser.add_argument('--middle', action = 'store_true', 
                        help = 'use middle point of each interval as the sentinel')
    parser.add_argument('--outer', action = 'store_true',
                        help = 'use outer point of each interval as the sentinel')
    parser.add_argument('--collapse', action = 'store_true', 
                        help = 'collapse initial interval to middle')
    # the collapsecenter is obsolete, instead use --collapse --outer
    # parser.add_argument('--collapsecenter', action = 'store_true',
    # help = 'collapse initial interval to middle but center one interval "on" the middle point')
    
    # controls flanking length
    # parser.add_argument('-f', '--flank', default=10000, type=int, 
    # help = 'length of flanking to plot, default 10kb')
    # --flank is also obsolete, it is replaced by --flankstep X --flanknumber Y
    parser.add_argument('--flanktoneighbor', action = 'store_true',
                        help = 'length of flanking is dependent on the nearest record')
    parser.add_argument('--flankbygene', action = 'store_true', # previously called varyflank
                        help = 'allow flanking region to vary according to the gene length')
    parser.add_argument('-f', '--flankstep', type = int, default=100, 
                        help = 'plot each X bases for flanking sequences, by default false (-1), this overrides -f')
    parser.add_argument('-m', '--flanknumber', type = int, default=30, 
                        help = 'number of points to sample in the flanking region')

    parser.add_argument('--expansion', type=int, default=0,
                        help = 'number of bases to expand in the two directions')
    # I don't think you need to specify window size, just do it by --middle and awk-expand
    # parser.add_argument('-w', '--windowsize', type = int, default=None,
    # help='Window Size, default to step size (non-overlapping windows)')

    # controls internal sampling
    parser.add_argument('-n', '--numinternal', type = int, default=30, 
                        help = 'number of points to sample in the genic/internal region, --middle ignores this')
    # others
    parser.add_argument('--fold', action = 'store_true', 
                        help = 'use the same index for intervals from two sides the target, usually used when strand is irrelevant')
    parser.add_argument('-s', '--strand', type=int, default=None, 
                        help = 'the field which contains strand information, if None then ignore strand')

    parser.set_defaults(func=main)
    args = parser.parse_args()
    try:
        args.func(args)
    except IOError:
        exit
