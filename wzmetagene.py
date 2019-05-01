#!/usr/bin/env python
"""
Generate intervals for meta-gene plot.
Usage:
python ~/wzlib/pyutils/wzmetagene.py input_bed | les
"""

import argparse
import numpy as np

def sample_forward(args, chrm, beg, end, fields, index_func):

    """ sample interval in the direction of bigger coordinates """

    _window_end = end + args.step
    for i in xrange(args.numflank):
        _window_beg = _window_end - args.windowsize
        if args.middle:
            window_mid = int((_window_beg + _window_end)/2.0)
            window_beg = window_mid-1
            window_end = window_mid
        else:
            window_beg = int(_window_beg)
            window_end = int(_window_end)
            _window_end += args.step 

        index = index_func(i)
        if index < 0:
            reg = '(%d)-(%d)' % (end - window_end, end - window_beg)
        else:
            reg = '%d-%d' % (window_beg - end, window_end - end)

        # in collapsed case, you can have index 0
        if index >= 0:
            area = 1
        else:
            area = -1
            
        if window_beg > 0 and window_end > window_beg:
            print('%s\t%d\t%d\t%d\t%s\t%d\t%s' %
                  (chrm, window_beg, window_end,
                   index, reg, area, '\t'.join(fields)))

def sample_backward(args, chrm, beg, end, fields, index_func):

    """ sample intervals in the direction of smaller coordinates """

    _window_beg = beg - args.step
    for i in xrange(args.numflank):
        _window_end = _window_beg + args.windowsize
        if args.middle:
            window_mid = int((_window_beg + _window_end)/2.0)
            window_beg = window_mid-1
            window_end = window_mid
        else:
            window_beg = int(_window_beg)
            window_end = int(_window_end)
            _window_beg -= args.step

        index = index_func(i)
        if index < 0:
            reg = '(%d)-(%d)' % (window_beg - beg, window_end - beg)
        else:
            reg = '%d-%d' % (beg - window_end, beg - window_beg)

        if index >= 0:
            area = 1
        else:
            area = -1
            
        if window_beg > 0 and window_end > window_beg:
            print('%s\t%d\t%d\t%d\t%s\t%d\t%s' %
                  (chrm, window_beg, window_end,
                   index, reg, area, '\t'.join(fields)))

def sample_internal(args, chrm, beg, end, fields, index_func):
    
    """ sample intervals internal of the input interval """

    sentinels = list(np.linspace(beg, end, args.numinternal+1))
    for i in xrange(len(sentinels)-1):
        window_beg = int(sentinels[i])
        window_end = int(sentinels[i+1])
        if args.middle:
            window_mid = int((sentinels[i] + sentinels[i+1])/2)
            window_beg = window_mid-1
            window_end = window_mid

        index = index_func(i)
        # 0 for internal
        if window_beg > 0 and window_end > window_beg:
            print('%s\t%d\t%d\t%d\t%d-%d%%\t0\t%s' %
                  (chrm, window_beg, window_end, index, 
                   float(index)/args.numinternal*100,
                   float(index+1)/args.numinternal*100, '\t'.join(fields)))

def main(args):

    if args.step > 0:
        args.flank = args.step * args.numflank
    elif args.numflank > 0:
        args.step = int(args.flank / args.numflank)

    if args.collapse:
        args.numinternal = 0

    if args.collapsecenter:
        args.numinternal = 1

    if args.windowsize is not None:
        args.windowsize = args.windowsize
    else:
        args.windowsize = args.step

    for line in args.table:
        fields = line.strip('\n').split('\t')
        chrm = fields[0]
        beg = int(fields[1])
        end = int(fields[2])

        if args.collapse:
            mid = (beg + end + 1) / 2
            beg = mid
            end = mid

        if args.collapsecenter:
            mid = (beg + end + 1) / 2
            beg = mid - args.windowsize / 2
            end = mid + args.windowsize / 2

        if args.strand is None:
            strand = '+'
        else:
            strand = fields[args.strand-1]

        if strand == '+':
            if args.fold:
                sample_backward(args, chrm, beg, end, fields, lambda i: -i-1)
                sample_internal(args, chrm, beg, end, fields, lambda i: min(i, args.numinternal-1-i))
                sample_forward(args, chrm, beg, end, fields, lambda i: -i-1)
            else:
                sample_backward(args, chrm, beg, end, fields, lambda i: -i-1)
                sample_internal(args, chrm, beg, end, fields, lambda i: i)
                sample_forward(args, chrm, beg, end, fields, lambda i: args.numinternal+i)
        else:
            sample_forward(args, chrm, beg, end, fields, lambda i: -i-1)
            sample_internal(args, chrm, beg, end, fields, lambda i: args.numinternal-1-i)
            sample_backward(args, chrm, beg, end, fields, lambda i: args.numinternal+i)

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate meta gene')
    parser.add_argument('table', help="Input bed file", type = argparse.FileType('r'), default='-')
    # report middle point for each sampled interval
    parser.add_argument('--middle', action = 'store_true', 
                        help = 'use middle point of each interval as the sentinel')
    parser.add_argument('--collapse', action = 'store_true', 
                        help = 'collapse initial interval to middle')
    parser.add_argument('--collapsecenter', action = 'store_true',
                        help = 'collapse initial interval to middle but center one interval "on" the middle point')
    # controls flanking length
    parser.add_argument('-f', '--flank', default=10000, type=int, 
                        help = 'length of flanking to plot, default 10kb')
    parser.add_argument('-F', '--step', type = int, default=-1, 
                        help = 'plot each X bases for flanking sequences, by default false (-1), this overrides -f')
    parser.add_argument('-m', '--numflank', type = int, default=30, 
                        help = 'number of points to sample in the flanking region')
    parser.add_argument('-w', '--windowsize', type = int, default=None, 
                        help='Window Size, default to step size (non-overlapping windows)')
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
