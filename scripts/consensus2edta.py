import csv, os, sys, pandas, numpy, re
try:
	from io import StringIO
except:
	from StringIO import StringIO

isomass = 1.0033548378

def read_consensus(fn):
    rd = csv.reader(open(fn), delimiter='\t')
    consHeader = []
    peptHeader = []
    runsName = []

    for row in rd:
        if row[0] == '#CONSENSUS':
            consHeader = row
        elif row[0] == '#UNASSIGNEDPEPTIDE':
            peptHeader = row
        elif row[0] == 'MAP':
            runsName.append(row[2])

    s = StringIO()
    with open(fn) as fh:
        for line in fh:
            if line.startswith("CONSENSUS"):
                s.write(line)
    s.seek(0)
    cons = pandas.read_csv(s, sep='\t', header=None, names=consHeader)
    
    coPeps = []
    with open(fn) as fh:
        for line in fh:
            if line.startswith("CONSENSUS"):
                coPeps.append('')
            elif line.startswith('PEPTIDE') and coPeps[-1] == '':
                coPeps[-1] = line.split("\t")[5]
    
    cons['peptide_0'] = coPeps
    
    s = StringIO()
    with open(fn) as fh:
        for line in fh:
            if line.startswith("UNASSIGNEDPEPTIDE"):
                s.write(line)
    s.seek(0)
    uaPeps = pandas.read_csv(s, sep='\t', header=None, names=peptHeader)
    uaPeps = uaPeps.groupby(['sequence', 'charge']).mean()
    return cons, uaPeps, runsName



if __name__ == '__main__':
    consensus = sys.argv[1]
    decoy_rt = float(sys.argv[2])
    decoy_ppm = float(sys.argv[3])
    out_dir = os.path.abspath(sys.argv[4])

    outpath = os.path.join(out_dir,  os.path.basename(os.path.splitext(consensus)[0]) + ".edta" )
    print(outpath)

    cons, uapep, runsName = read_consensus(consensus)
    with open(outpath, 'w') as fh:
        fh.write("\t".join(['RT', 'MZ', 'Int', 'Charge', 'rank']) + '\n')
        rank = 0
        for i, x in cons.iterrows():
            rank += 1
            s = x.peptide_0
            s = re.sub('C\(Carbamidomethyl\)', 'C', s)
            z = x.charge_cf
            rt = x.rt_cf
            mz = x.mz_cf
            ion = x.intensity_cf
            

            # Target 
            fh.write( "%.3f\t%f\t%f\t%s\t%.1f\n" % (rt, mz, ion, z, rank + 0.0) )
            mz1 = isomass / z + mz
            fh.write( "%.3f\t%f\t%f\t%s\t%.1f\n" % (rt, mz1, ion, z, rank+ 0.1) )
            
            # Decoys 
            mz += numpy.random.normal(decoy_ppm * mz * 1e-6 , mz * 1e-6)
            rt += numpy.random.normal(decoy_rt, 5)

            fh.write( "%.3f\t%f\t%f\t%s\t%.1f\n" % (rt, mz, ion, z, rank+ 0.2) )
            mz1 = isomass / z + mz
            fh.write( "%.3f\t%f\t%f\t%s\t%.1f\n" % (rt, mz1, ion, z, rank+ 0.3) )


        for i, x in uapep.iterrows():
            rank += 1
            s, z = i
            rt = x.rt
            mz = x.mz

            fh.write("%.3f\t%f\t0\t%s\t%.1f\n" % (rt, mz, z, rank + 0.0))
            mz = isomass / z + mz
            fh.write("%.3f\t%f\t0\t%s\t%.1f\n" % (rt, mz, z, rank + 0.1))

            mz += numpy.random.normal(decoy_ppm * mz * 1e-6 , mz * 1e-6)
            rt += numpy.random.normal(decoy_rt, 5)

            fh.write( "%.3f\t%f\t0\t%s\t%.1f\n" % (rt, mz, z, rank + 0.2) )
            mz1 = isomass / z + mz
            fh.write( "%.3f\t%f\t0\t%s\t%.1f\n" % (rt, mz1, z, rank + 0.3) )


