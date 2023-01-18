#!/usr/bin/env python3

import argparse
import cvlrcommon
import matplotlib.pylab as plt
import numpy as np
import re
from scipy.special import betaln
import sys


def h(a1, b1, a2, b2):
    tmp = betaln(a1+a2, b1+b2) - betaln(a1,b1) - betaln(a2,b2)
    return(np.exp(tmp))

def g(a1, b1, a2, b2):
    if ( a1 == a2 and b1 == b2):
        return 0.5
    elif  (a1>a2):
        return g(a1 - 1, b1, a2, b2)  + h (a1 - 1, b1, a2,  b2) / ( a1 - 1 )
    elif (a2>a1):
        return g(a1,  b1, (a2 - 1), b2)  - h (a1, b1, (a2 - 1), b2 ) / ( a2 - 1 )
    elif (b1>b2):
        return g ( a1,  (b1 - 1),  a2, b2 ) -  h( a1,  (b1 - 1),  a2,  b2) / ( b1 - 1 ) 
    elif (b2>b1):
        return g(  a1, b1, a2, ( b2 - 1 ))  + h (a1, b1, a2, ( b2 - 1 ))  / ( b2 - 1 )
    else:
        raise ValueError

def betadiff(n01, n11, n02, n12):
    return g((n11+1), (n01+1), (n12+1), (n02+1))
    
def n0_n1(n, d, k, dstate, drnames, cl):
    """
    returns 2 arrays n0, n1 (k,d)
    which store the number of 0s (or 1s) for cluster
    k at position d.
    This is used to compute the mean methylation per cluster
    """
    n0 = np.full((k, d), fill_value = 0 )
    n1 = np.full((k, d), fill_value = 0 )
    notfound = set()
    found = 0
    ## the computation below uses two dictionaries
    ## dstate and cl.
    ## dstate is built by parsing the matrix
    ## cl is built by parsing the clusters (or the output of whatshap).
    ## There might be keys in dstate which do not appear in cl,
    ## and vice versa.
    for (i, j) in dstate.keys():
        try:
            if ( 0 == dstate[(i, j)] ):
                n0[cl[drnames[i]], j] += 1
            if ( 1 == dstate[(i, j)] ):
                n1[ cl[drnames[i]], j ] += 1
        except KeyError:
            notfound.add(i)
            continue
    return (n0, n1, notfound)

def parse_info(args):
    clusterfn = args.clusterfn
    clusterf = open(clusterfn, 'r')
    k = -1; mu=[]
    for line in clusterf:
        line = line.strip()
        if (re.match("#@K:", line)):
                line=re.sub("#@K:","",line)
                k = int(line)
                print(f"K={k}")
        if (re.match("#@GPOS:", line)):
                print("GPOS found")
                line = re.sub("#@GPOS:","",line)
                gpos = [int(e) for e in line.split("\t")]
        if (re.match("#@PI:", line)):
                line = re.sub("#@PI:","",line)
                pi = [float(e) for e in line.split("\t")]
        for i in range(k):
                mstr = f"#@MU{i}:"
                if (re.match(mstr, line)):
                        print(f"MU{i} found")
                        line=re.sub(mstr,"",line)
                        mu.append([float(e) for e in line.split("\t")])
    clusterf.close()
    return (k, np.array(gpos), np.array(pi), np.array(mu))

def parse_haplo(haplofn):
    f = open(haplofn, 'r')
    pset = set()
    # cl is a dictionary rname -> cluster id
    cl={} 
    nr = np.full(2, fill_value=0)
    for line in f:
        if ('#' == line[0]): continue
        line = line.strip()
        fields = line.split("\t")
        rname = fields[0]
        haplo = fields[1]
        pset.add(fields[2])
        if ( "H1" == haplo ):
            cl[rname]=0
            nr[0]+=1
        if ( "H2" == haplo ):
            cl[rname]=1
            nr[1]+=1
    f.close()
    return (nr, cl, pset)

def compute_stats(d, k, n0, n1):
    """
    return:
    avgmeth (pos, k) methylation at position pos for cluster k
    meth (pos) methylation at position pos
    diff pos, k(k-1)/2 absolute difference in methylation between clusters at position pos
    postprob pos,k(k-1)/2 array i,j posterion probability that methylation of cluster i > meth at cluster j at position pos 
    """
    avgmeth = np.full((d,k), fill_value=0.0, dtype=np.float64)
    meth = np.full(d, fill_value=0, dtype=np.float64)
    diff = np.full((d,k*(k-1)//2), fill_value=0, dtype=np.float64)
    postprob = np.full((d,k*(k-1)//2), fill_value=0, dtype=np.float64)
    for j in range(d):
        for l in range(k):
            avgmeth[j,l] = n1[l, j] / (n0[l, j] + n1[l, j])
        meth[j] = np.sum(n1[:,j])/(np.sum(n0[:,j])+np.sum(n1[:,j]))
        didx=0
        for l1 in range(k):
            for l2 in range(l1+1,k):
                diff[j, didx] = np.abs(avgmeth[j,l1] - avgmeth[j,l2])
                postprob[j, didx] = betadiff( n0[l1, j], n1[l1,j],n0[l2,j],n1[l2,j]  )
                didx += 1
    return (avgmeth, meth, diff, postprob)

def median_abs_diff_notna(diff):
    """
    return array with median abs diff for each pair of clusters
    """
    n, npairs = diff.shape
    k = int(1+np.sqrt(1+8*npairs))//2
    madnna = np.full((k,k), fill_value = -1, dtype=float)
    diffidx=0
    for l1 in range(k):
        for l2 in range(l1+1,k):
            notnaidx = np.isfinite(diff[:,diffidx])
            madnna[l1,l2] = np.median(diff[notnaidx, diffidx])
            diffidx+=1
    return madnna
    
def print_stats(n0, n1, avgmeth, meth, diff,postprob, dgpos, dmrt):
    """
    print statistics and optionally DMRs.
    With the option --dmrs the output has an extra column for each cluster pair (after the posterior probability)
    This will contain either '-' or 'DMR_c1_c2_\d+' indicating whether
    the corresponding position belongs to some DMR between clusters c1 and c2.
    The DMR coordinates are interspersed with the normal output with the
    format:
    #@DMR_c1_c2_\d+:<DMR start, stop, median abs diff>
    """
    d,k = avgmeth.shape
    madnna = median_abs_diff_notna(diff)
    for l1 in range(k):
        for l2 in range(l1+1,k):
            print(f"#@MEDIAN_ABS_DIFF_NOTNA_{l1}_{l2}:{madnna[l1,l2]:.3f}")
    if ( -1 == dmrt ):
        dmrt = 2
    ndmr = np.full((k,k), fill_value = 0)
    indmr = np.full((k,k), fill_value=False)
    for j in range(d):
        print(f"{dgpos[j]}", end="\t")
        for l in range( k ):
            print(f"{n0[l, j]}\t{n1[l, j]}\t{avgmeth[j,l]:.3f}", end="\t")
        sum0 = np.sum(n0[:,j])
        sum1 = np.sum(n1[:,j])
        depth = sum0 + sum1
        print(f"{sum0}\t{sum1}\t{meth[j]:.3f}\t{depth}", end="")
        l = 0
        ldmrs = []
        for l1 in range(k):
            for l2 in range(l1+1,k):
                print(f"\t{diff[j,l]:.3f}", end="")
                print(f"\t{postprob[j,l]:.3g}", end="")
                ## DMRs are decided looking at the folded posterior probability
                fpp = postprob[j,l] if ( postprob[j,l] > 0.5 ) else 1.0 - postprob[j,l] 
                if ( fpp > dmrt ):
                    if ( not indmr[l1,l2] ):
                        ndmr[l1,l2] +=1
                    msg = f"\tDMR_{l1}_{l2}_{ndmr[l1,l2]}" 
                    indmr[l1,l2] = True
                else:
                    if (indmr[l1,l2]):
                        ldmrs.append((l1,l2, ndmr[l1,l2]))
                    indmr[l1,l2] = False
                    msg="\t-"
                if (dmrt != 2):
                    print(msg, end="")
            l+=1
        print()
        ## here print the stored dmrs        
        for e in ldmrs:
            print(f"#@{e}")
        
def covmatrix(mu, pi):
    k = mu.shape[0]
    d = mu.shape[1]
    pi = pi.reshape((k, 1))
    cov = np.full((d, d), fill_value=0.)
    expx = np.matmul(pi.T, mu)
    for l in range(k):
        sigma = np.diag(mu[l, :]*(1-mu[l, :]))
        cov = cov + pi[l] * sigma
        rl = mu[l, :].reshape((1, d))
        cov = cov + np.matmul(rl.T, rl)
    cov = cov - np.matmul(expx.T, expx)
    return cov

def random_clusters(drnames):
    k = 2
    nr = np.full(k, fill_value=0)
    cl = {}  # rname ==> clusterid
    for rname in drnames.values():
        clidx = np.random.choice(k)
        nr[clidx] += 1
        cl[rname] = clidx
    ##print(cl)
    return (nr,cl)

#####################
# commands          #
#####################
        
def haplo_or_cluster(args):
    matrixfn  = args.matrixfn
    dmrt = float(args.dmrt)
    drnames, dgpos, dstate, maxridx, maxgposidx = cvlrcommon.gmatrix_of_file(matrixfn)
    n = maxridx + 1
    d = maxgposidx + 1
    cmd = ""
    try:
        haplofn = args.haplofn
        (nr, cl, pset) = parse_haplo(haplofn)
        cmd = "haplo"
    except AttributeError:
        clusterfn = args.clusterfn
        (nr, cl) = cvlrcommon.parse_clusters(clusterfn)
        cmd = "cluster"
    k = len(nr)
    print(f"#@N:{n}")
    print(f"#@K:{k}")
    for i in range(k):
        print(f"#@POP{i}:{nr[i]}")
    ( n0, n1 , notfound ) = n0_n1(n, d, k, dstate, drnames, cl)
    print(f"#@NOTFOUND:{len(notfound)}")
    print(f"#@CMD:{cmd}")
    if ( "haplo" == cmd ):
        print(f"#@HAPLOFN:{haplofn}")
    elif ( "cluster" == cmd):
        print(f"#@CLUSTERFN:{clusterfn}")
    else:
        print("invalid command :{cmd}", file=sys.stderr)
        sys.exit(1)
    print(f"#@MATRIXFN:{matrixfn}")
    ( avg, meth, diff, postprob ) = compute_stats(d, k, n0, n1)
    print_stats(n0, n1,avg,meth, diff, postprob,dgpos,dmrt)

def mean(args):
    pass

def cov(args):
    (k, gpos, pi, mu) = parse_info(args)
    covm = covmatrix(mu, pi)
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes()
    im = ax.imshow(covm)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel( "covariance", rotation=-90, va="bottom" )
    lstep = 300
    labcol= [f"{gpos[e]}" for e in range(lstep,len(gpos),lstep)]
    ax.set_xticks(np.arange(lstep, covm.shape[1], lstep), labels=labcol)
    ax.tick_params(axis='x', labelright=True, labelrotation=-30)
    ax.set_yticks(np.arange(lstep, covm.shape[1], lstep), labels=labcol)
    if (args.outputfn != None):
        plt.savefig(args.outputfn)
    else:
        plt.show()

def randomize(args):
    matrixfn = args.matrixfn
    nsamples = int(args.nsamples)
    drnames, dgpos, dstate, maxridx, maxgposidx = cvlrcommon.gmatrix_of_file(matrixfn)
    n = maxridx + 1
    d = maxgposidx + 1
    k = 2
    for s in range(nsamples):
        print(f"#@N:{n}")
        print(f"#@K:{k}")
        (nr, cl) = random_clusters(drnames)
        for i in range(k):
            print(f"#@POP{i}:{nr[i]}")
        ( n0, n1 , notfound ) = n0_n1(n, d, k, dstate, drnames, cl)
        print(f"#@NOTFOUND:{notfound}")
        ( avg, meth, diff, postprob ) = compute_stats(d, k, n0, n1)
        print_stats(n0, n1,avg,meth, diff, postprob, dgpos)

def pval(args):
    clusterfn = args.clusterfn
    matrixfn = args.matrixfn
    nsamples = int(args.nsamples)
    drnames, dgpos, dstate, maxridx, maxgposidx = cvlrcommon.gmatrix_of_file(matrixfn)
    n = maxridx + 1
    d = maxgposidx + 1
    k = 2
    ## true value
    (nr, cl) = cvlrcommon.parse_clusters(clusterfn)
    ( n0, n1 , notfound ) = n0_n1(n, d, k, dstate, drnames, cl)
    ( avg, meth, diff ) = compute_stats(d, k, n0, n1)
    madnna_true = median_abs_diff_notna(diff)[0,1]
    print(f"#@MATRIXFN:{matrixfn}")
    print(f"#@MEDIAN_ABS_DIFF_NOTNA:{madnna_true:.3f}")
    samples = np.full((nsamples,), fill_value=-1, dtype=float)
    for s in range(nsamples):
        (nr, cl) = random_clusters(drnames)
        ( n0, n1 , notfound ) = n0_n1(n, d, k, dstate, drnames, cl)
        ( avg, meth, diff ) = compute_stats(d, k, n0, n1)
        madnna = median_abs_diff_notna(diff)
        samples[s] = madnna[0,1]
        print(f"sample[{s}]:{samples[s]}", file=sys.stderr)
    pval = np.sum(samples > madnna_true) / nsamples
    print(f"#@PVAL_MADNNA:{pval:.3g}")
        
    
parser     = argparse.ArgumentParser("cvlr-stats")
subparsers = parser.add_subparsers(help="sub-command help")
# haplo
parser_haplo = subparsers.add_parser( "haplo", help="compute methylation on whatshap haplotypes" )
parser_haplo.add_argument( "haplofn",  help="haplo file" )
parser_haplo.add_argument( "matrixfn",  help="matrix file" )
parser_haplo.set_defaults( func = haplo_or_cluster )
# cluster
parser_cluster = subparsers.add_parser("cluster", help="compute methylation on clusters")
parser_cluster.add_argument("clusterfn",  help = "cluster file (output of cvlr-cluster)")
parser_cluster.add_argument("matrixfn",  help="matrix file")
parser_cluster.add_argument("--dmrs", help="threshold for optional dmr computation",  dest = "dmrt", default=-1)
parser_cluster.set_defaults(func = haplo_or_cluster)
# mean
parser_mean = subparsers.add_parser("mean", help="extracts the mean per cluster as inferred by cvlr-cluster")
parser_mean.add_argument("clusterfn",  help = "cluster file (output of cvlr-cluster)")
parser_mean.set_defaults(func=mean)
# cov
parser_cov = subparsers.add_parser("cov", help="computes the covariance matrix from the output of cvlr-cluster")
parser_cov.add_argument("clusterfn",  help = "cluster file (output of cvlr-cluster)")
parser_cov.add_argument("--outputfn", help = "don't show plot, save in file")
parser_cov.set_defaults(func = cov)
# randomize
parser_randomize = subparsers.add_parser("randomize", help="computes stats over randomized (binary) clusters")
parser_randomize.add_argument("matrixfn",  help="matrix file")
parser_randomize.add_argument("nsamples",  help="number of samples")
parser_randomize.set_defaults(func = randomize)
# pval
parser_pval = subparsers.add_parser("pval", help="computes pval of median abs diff (k=2)")
parser_pval.add_argument("clusterfn",  help="cluster file")
parser_pval.add_argument("matrixfn",  help="matrix file")
parser_pval.add_argument("nsamples",  help="number of samples")
parser_pval.set_defaults(func = pval)

if (__name__ == '__main__'):
    args = parser.parse_args()
    try:
        func = args.func
    except AttributeError:
        parser.parse_args(['-h'])
    func(args)
    
