# DeMix-Q. Peptide quantification
# bo.zhang@ki.se

from __future__ import division, print_function
import pandas
import numpy as np
import sys
import csv
import re
from six import StringIO
from sklearn.neighbors import KNeighborsRegressor
from scipy.stats import pearsonr


def read_consensus(fn):
    '''
    Read text table from consensusXML exported by OpenMS TextExporter
    '''
    cons_header = []
    pept_header = []
    runs_name = []

    # fetch the headers for consensus features, unassigned peptides and
    # experiments' names.
    for row in csv.reader(open(fn), delimiter='\t'):
        if row[0] == '#CONSENSUS':
            cons_header = row
        elif row[0] == '#UNASSIGNEDPEPTIDE':
            pept_header = row
        elif row[0] == 'MAP':
            runs_name.append(row[2])

    # read consensus features
    s = StringIO()
    with open(fn) as fh:
        for line in fh:
            if line.startswith("CONSENSUS"):
                s.write(line)
    s.seek(0)
    cons = pandas.read_csv(s, sep='\t', header=None, names=cons_header)
    co_peps = []
    with open(fn) as fh:
        for line in fh:
            if line.startswith("CONSENSUS"):
                co_peps.append('')
            elif line.startswith('PEPTIDE') and co_peps[-1] == '':
                # choose the first recorded peptide sequence as consensus
                # sequence
                co_peps[-1] = line.split("\t")[5]
    cons['peptide_0'] = co_peps

    # read uassigned peptides as consensus features
    s = StringIO()
    with open(fn) as fh:
        for line in fh:
            if line.startswith("UNASSIGNEDPEPTIDE"):
                s.write(line)
    s.seek(0)
    ua_peps = pandas.read_csv(s, sep='\t', header=None, names=pept_header)
    ua_peps = ua_peps.groupby(['sequence', 'charge']).mean()

    return cons, ua_peps, runs_name


def read_eic(fn):
    '''
    Read the detailed output file from EIC extraction in OpenMS.
    Two isotopic peaks (M and M+1) are extracted for each consensus feature and unassigned peptide
    Calculate geometric average of intensities from the two isotopic peaks, as well as diviations of RT and mass.
    '''

    with open(fn) as fh:
        eic_str = [StringIO() for _ in range(4)]
        sample_header = [i for i in next(fh).rstrip().split(',')
                         if i]  # sample row
        next(fh)  # empty row
        cols = next(fh).rstrip().split(',')  # quantity headers
        for ix in range(2, len(cols)):
            i = int((ix - 2) / 5)
            # rename columns according to sample names
            cols[ix] = '_'.join([sample_header[i], cols[ix]])
        ix = 0
        for line in fh:
            eic_str[ix % 4].write(line)
            ix += 1
        [sio.seek(0) for sio in eic_str]

    # obtain quantities from M and M+1 of target and decoy features

    # Monoisotopic
    eic = pandas.read_csv(eic_str[0], header=None, names=cols)
    # 13C isotope
    eic_iso = pandas.read_csv(eic_str[1], header=None, names=cols)

    # Monoisotopic (decoy)
    eic_decoy0 = pandas.read_csv(eic_str[2], header=None, names=cols)
    # 13C isotope (decoy)
    eic_decoy1 = pandas.read_csv(eic_str[3], header=None, names=cols)

    for samp in sample_header:
        int_ix = samp + '_intensity'
        rt_ix = samp + '_dRT'
        ppm_ix = samp + '_dppm'

        # values of target features
        eic[samp + '_int_1'] = (eic[int_ix] * eic_iso[int_ix]) ** 0.5
        eic[samp + '_dRT_1'] = eic[rt_ix] - eic_iso[rt_ix]
        eic[samp + '_ppm_1'] = eic[ppm_ix] - eic_iso[ppm_ix]

        # values of decoy features
        eic[samp + '_int_d1'] = (eic_decoy0[int_ix] *
                                 eic_decoy1[int_ix]) ** 0.5
        eic[samp + '_dRT_d0'] = eic_decoy0[rt_ix]
        eic[samp + '_dRT_d1'] = eic_decoy0[rt_ix] - eic_decoy1[rt_ix]
        eic[samp + '_ppm_d0'] = eic_decoy0[ppm_ix]
        eic[samp + '_ppm_d1'] = eic_decoy0[ppm_ix] - eic_decoy1[ppm_ix]

    return eic


def median_shift(dat, rcol, scol, knn=15, span=500):
    '''
    Predict median shift (fold-change) of ion intensity at a given RT, using KNN regression
    rcol: column of global medians
    scol: columm of extracted intensities of one sample
    knn:  number of k nearest neighbors
    span: number of features for median shift calculation.
    '''
    d = dat[(dat[scol] > 0) & (dat[rcol] > 0)].sort_values(by='rt_cf')
    x = d['rt_cf'].values
    # fold-change
    y = np.log2(d[scol] / d[rcol])
    X = np.array([np.mean(x[i:i + span])
                  for i in range(0, len(x), int(span / 10))])
    y = np.array([np.median(y[i:i + span])
                  for i in range(0, len(x), int(span / 10))])
    reg = KNeighborsRegressor(n_neighbors=knn)
    reg.fit(np.matrix(X).T, np.matrix(y).T)
    return reg, X, y


def main():
    import argparse
    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('-consensus', default='consensus.csv',
                       help='The consensus feature map')

    apars.add_argument('-eic', default='eic.csv',
                       help='The extracted ion-chromatography')

    apars.add_argument('-n_samples', default='4', help='number of samples')
    apars.add_argument('-n_repeats', default='3',
                       help='number of replicate experiments per sample')
    apars.add_argument('-fdr', default='0.05',
                       help='Quality cutoff by estimating FDR from decoy extractions')

    apars.add_argument('-knn_k', default='5',
                       help='K for KNN correction, set 0 to disable correction')
    apars.add_argument('-knn_step', default='500',
                       help='Number of sequential features used for estimating local medians')

    apars.add_argument('-out', default='peptides.DeMixQ.csv',
                       help='Output table in csv format')

    # -----------------------
    args = apars.parse_args()
    consensus_txt = args.consensus
    eic_txt = args.eic
    num_samples = int(args.n_samples)
    num_replica = int(args.n_repeats)
    fdr = float(args.fdr)
    knn_k = int(args.knn_k)
    knn_step = int(args.knn_step)
    # -----------------------
    eic = read_eic(eic_txt)
    cons, uapep, ic0 = read_consensus(consensus_txt)
    icols = [i for i in cons.columns if i.startswith('intensity_')]
    cons = cons[icols + ['quality_cf', 'peptide_0',
                         'charge_cf', 'rt_cf', 'mz_cf']]
    cons.rename(columns=dict(zip(icols[1:], ic0)), inplace=True)

    # please pardon me for the confusing variable namings here.
    # ic1x => feature Intensity Column 1 :XIC (geometric average of M and M+1)
    ic1x = [i for i in eic.columns if i.endswith('int_1')]

    # column names of consensus (targets)
    crt0 = [i for i in eic.columns if i.endswith('dRT')]
    crt1 = [i for i in eic.columns if i.endswith('dRT_1')]
    cppm0 = [i for i in eic.columns if i.endswith('ppm')]
    cppm1 = [i for i in eic.columns if i.endswith('ppm_1')]

    # column names of decoys
    dic1x = [i for i in eic.columns if i.endswith('int_d1')]
    drt0 = [i for i in eic.columns if i.endswith('dRT_d0')]
    drt1 = [i for i in eic.columns if i.endswith('dRT_d1')]
    dppm0 = [i for i in eic.columns if i.endswith('ppm_d0')]
    dppm1 = [i for i in eic.columns if i.endswith('ppm_d1')]

    # combine Feature and XIC into dataframe_X
    dx = pandas.concat([eic, cons], axis=1)
    dx['peptide_0'] = [i for i in cons.peptide_0] + [i[0]
                                                     for i in uapep.index.tolist()]
    dx['charge_cf'] = [i for i in cons.charge_cf] + [i[1]
                                                     for i in uapep.index.tolist()]
    dx['mz_cf'] = [i for i in cons.mz_cf] + uapep.mz.tolist()
    dx['rt_cf'] = [i for i in cons.rt_cf] + uapep.rt.tolist()
    dx = dx[dx.rt_cf > 0]

    dx['peptide'] = [re.sub('C\(Carbamidomethyl\)', 'C', str(i))
                     for i in dx.peptide_0]
    dx['baseseq'] = [re.sub('\(.+?\)', '', str(i)) for i in dx.peptide_0]

    dx['mods'] = [sorted(re.findall('\(.+?\)', str(i))) for i in dx.peptide]
    dx['uniq'] = ["%s%d%s" % (x.baseseq, x.charge_cf, ''.join(x.mods))
                  for i, x in dx.iterrows()]

    # cross-run quantifications from feature-based linking.
    dx['f_overlap'] = [np.count_nonzero(
        np.nan_to_num(i)) for i in dx[ic0].values]
    # cross-run quantifications from ion-based linking.
    dx['e_overlap'] = [np.count_nonzero(i) for i in dx[ic1x].values]
    # median intensity in each individual run
    dx['medianEIC'] = dx[ic1x].median(axis=1)

    print("number of features:\t", len(dx))
    print("number of runs:\t", len(ic1x))

    # XIC extraction rate for consensus features
    e = np.array(dx[dx.intensity_cf.isnull() == False].e_overlap).tolist()
    print(e.count(len(ic1x)) * 100 / len(e), "% of ", len(e),
          " consensus features are extracted across all runs.")

    # Reference set: consensus features with no missing values in both feature
    # linking and XIC extraction
    refdx = dx[(dx.f_overlap == len(ic0)) & (dx.e_overlap == len(ic0))]

    tcol = np.concatenate([crt0, crt1, cppm0, cppm1])
    dcol = np.concatenate([drt0, drt1, dppm0, dppm1])

    # unit-less transformation
    refz = (refdx[tcol] - refdx[tcol].mean()) / refdx[tcol].std()
    decoyz = (refdx[dcol] - refdx[tcol].mean().values) / \
        refdx[tcol].std().values
    testz = (dx[tcol] - refdx[tcol].mean()) / refdx[tcol].std()

    for i in np.array(ic1x).reshape(num_samples, num_replica):
        cv = refdx[i].std(axis=1) / refdx[i].mean(axis=1)
        cv = np.sqrt(cv)
        for run_id in i:
            refz[run_id + '_cv'] = cv

    for i in np.array(dic1x).reshape(num_samples, num_replica):
        cv = refdx[i].std(axis=1) / refdx[i].mean(axis=1)
        cv = np.sqrt(cv)
        for run_id in i:
            decoyz[run_id + '_cv'] = cv

    for i in np.array(ic1x).reshape(num_samples, num_replica):
        cv = dx[i].std(axis=1) / dx[i].mean(axis=1)
        cv = np.sqrt(cv)
        for run_id in i:
            testz[run_id + '_cv'] = cv

    ccv = [i for i in refz.columns if i.endswith('cv')]
    dcv = [i for i in decoyz.columns if i.endswith('cv')]

    print('Average mass errors:\n', refdx[cppm0].mean())
    print('SD of mass errors:\n', refdx[cppm0].std())

    fil_cols = list(zip(ic1x, crt0, crt1, cppm0, cppm1, ccv))
    fil_cols_decoy = list(zip(dic1x, drt0, drt1, dppm0, dppm1, dcv))

    zscores = []
    for f in fil_cols:
        e = refdx[f[0]] > 0
        rz = refz[e]
        z = [np.sum(v * v) for v in rz[list(f[1:])].values]
        zscores = zscores + z

    dzscores = []
    for f in fil_cols_decoy:
        e = refdx[f[0]] > 0
        dd = decoyz[e]
        z = [np.sum(v * v) for v in dd[list(f[1:])].values]
        dzscores = dzscores + z

    score_list = sorted(zip(zscores + dzscores,
                            [0 for i in zscores] + [1 for i in zscores]),
                        key=lambda x: x[0])

    score_cutoff = 0
    hit_count = 0
    decoy_count = 0

    # set FDR threshold
    for s in score_list:
        hit_count += 1
        score_cutoff = -np.log(s[0])
        if s[1] > 0:
            decoy_count += 1
        if decoy_count / hit_count >= fdr:
            print("\nFDR(%.2f) cutoff score: %.2f" % (fdr, score_cutoff))
            break

    for f in fil_cols:
        score = np.array([np.sum(v * v) for v in testz[list(f[1:])].values])
        score = -np.log(score)
        score[np.isnan(score)] = -np.inf
        fraction = sum(score > score_cutoff) * 100. / sum(np.isfinite(score))
        print(f[0], '%.2f' % fraction, '%',
              'Feature-EIC pairs passed the score cutoff.')
        dx.loc[score <= score_cutoff, f[0]] = 0  # XIC filter
        # print((score > score_cutoff).sum() / len(score))

    dx['medianEIC'] = dx[dx[ic1x] > 0][ic1x].median(axis=1)
    dx['e_overlap'] = [np.count_nonzero(i) for i in dx[ic1x].values]

    print("Features before filtering:\t", len(dx))

    # quantified at least in one run in each sample
    ix = np.min([np.mean(dx[i].values, axis=1) for i in np.array(
        ic1x).reshape(num_samples, num_replica)], axis=0) > 0

    # or quantified at least in half of samples
    ix = dx.e_overlap > 0.5 * num_samples * num_replica

    # dx: missing value filtered dataframe_X
    dx = dx[ix]

    print("Features after filtering:\t", len(dx))
    print("Unique peptide sequences:\t", len(dx.baseseq.unique()), '\n')

    icols = list(zip(ic0, ic1x))
    print('Experiment', 'Features', 'Feature-EIC-correlation', sep='\t')
    for a in icols[:]:
        do = dx[pandas.notnull(dx[a[0]])]
        do = do[do[a[1]] > 0]
        xy = do[list(a)].apply(np.log2).values
        coeff = pearsonr(xy[:, 0], xy[:, 1])[0]
        print(a[1], len(xy), "%.2f" % coeff, sep='\t')
        if coeff < 0.7:
            print(a, " has weak correlation! Check the RT alignment result!")
            # icols.remove(a)

    # log2 transform of selected intensity colums from the dataframe_X
    data = dx[ic0 + ic1x].apply(np.log2)

    # impute remaining missing values as the minimum detectable intensity
    # data[data == -np.inf] = data[data != -np.inf][ic1x].min().min()

    # or impute missing values as zero intensity
    data[data == -np.inf] = 0

    data['rt'] = dx['rt_cf']
    data['mass'] = dx['mz_cf'] * dx['charge_cf']
    data['charge'] = dx['charge_cf']
    data['peptide'] = dx.peptide
    data['baseseq'] = dx.baseseq
    data['uniq'] = dx.uniq

    for a, b in icols:
        dn = data[(data[a].notnull()) & (data[b] > 0)]
        mx = np.matrix(dn[[a, b, 'mass']])

        # KNN regression (k=5 by default) for predicting feature abundance
        # based on ion intensity and precosor mass
        regr = KNeighborsRegressor().fit(mx[:, 1:], mx[:, 0])
        a_ = regr.predict(np.matrix(data[[b, 'mass']]))[:, 0]
        a_[data[b].values == 0] = 0  # keep missing values from extraction
        data[a + '_'] = a_

    ic2 = [a + '_' for a in ic0]
    tmp = dx[ic0].apply(np.log2)
    tmp.values[np.isnan(tmp.values)] = data[ic2].values[np.isnan(tmp.values)]
    data[ic0] = tmp.apply(np.exp2)  # converting back to linear scale

    fraction = (data[ic0].values == 1).sum() * 100 / data[ic0].values.size
    print('\nFraction of remaining missing values: %.2f' % fraction, r"%")

    fraction = len(data[data.peptide != '']) * 100. / data[ic0].__len__()
    print("\n%.2f" % fraction, r'%',
          "features have been assigned with peptide sequences.")

    # median-shift correction
    if knn_k:
        print('Correction of feature abundances\nK=%d' % knn_k,
              '\nStep=%d' % knn_step, 'features\n')
        for a, b in icols:
            reg = median_shift(dx, 'medianEIC', b, knn=knn_k, span=knn_step)[0]
            data[a] = np.log2(data[a]) - reg.predict(np.matrix(data.rt).T).T[0]
        data[ic0] = 2 ** data[ic0]

    # ========
    # Report peptide quantification results
    # keep identified features
    pepdata = data[data.peptide != ''].drop_duplicates(subset=['uniq'])

    # remove modified peptides?
    # pepdata = pepdata[[len(i) == 0 for i in pepdata.mods]]

    mx = pepdata.groupby(['baseseq'])[ic0].sum().to_csv(args.out)

    print('Finished.\nOutput csv: ', args.out)

if __name__ == '__main__':
    main()
