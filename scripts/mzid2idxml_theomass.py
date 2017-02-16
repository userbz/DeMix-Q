from __future__ import print_function
import sys
import numpy
import os
import lxml.etree as et
from pyteomics import mzid
from pprint import pprint
# Bo Zhang 2015.01.12 : Use theoretical m/z
# Python3


def _cleanup(e):
    e.clear()
    while e.getprevious() is not None:
        del e.getparent()[0]


def _dbrefs(fn):
    ns = '{http://psidev.info/psi/pi/mzIdentML/1.1}'
    refDict = dict()
    protRef = dict()
    xml = et.iterparse(fn, events=('start', 'end'))
    for i, e in xml:
        if e.tag == ns + 'DBSequence' and i == 'end':
            protRef[e.get('id')] = e.get('accession')
            _cleanup(e)
        elif e.tag == ns + 'Peptide' and i == 'end':
            pepid = e.get('id')
            baseseq = next(e.iter(ns + 'PeptideSequence')).text
            peplist = [''] + list(baseseq)
            for mod in e.iter(ns + 'Modification'):
                pos = int(mod.get('location'))
                unimod = mod.getchildren()[0].get('name')
                peplist[pos] = peplist[pos] + '(%s)' % unimod
            pepseq = ''.join(peplist)
            refDict[pepid] = (baseseq, pepseq)
            _cleanup(e)
        elif e.tag == ns + 'PeptideEvidence' and i == 'end':
            flag = e.get('isDecoy') == 'false' and 1 or -1
            refDict[e.get('id')] = (flag, e.get('pre'), e.get('post'), protRef[e.get('dBSequence_ref')])
            _cleanup(e)
        elif e.tag == ns + 'SequenceCollection' and i == 'end':
            del xml
            return refDict


def mzid2pin(fn, max_charge=7):
    fh = open(fn + '.idxml', 'w')
    db = _dbrefs(fn)
    psmReader = mzid.read(fn)
    print( r'''<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="http://open-ms.sourceforge.net/XSL/IdXML.xsl" ?>
<IdXML version="1.2" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/SCHEMAS/IdXML_1_2.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <SearchParameters id="SP_0" db="" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="Trypsin/P" missed_cleavages="3" precursor_peak_tolerance="10" peak_mass_tolerance="20" ></SearchParameters>
    <IdentificationRun date="" search_engine="" search_engine_version="" search_parameters_ref="SP_0" > ''', file=fh)

    for psm in psmReader:
        if 'SpectrumIdentificationItem' not in psm:
            continue
        rt = psm['scan start time']
        for match in psm['SpectrumIdentificationItem']:
            if match['MS-GF:QValue'] > 0.01:
            # if match['MS-GF:PepQValue'] > 0.01:
                continue
            z = match['chargeState']
            if not 2 <= z <= max_charge:
                continue
            expMZ = match['experimentalMassToCharge']
            mz = match['calculatedMassToCharge']  # - (expMZ - featMZ)
            score = match['MS-GF:EValue']

            flag, enzN, enzC, prot = db[match['PeptideEvidenceRef'][0]['peptideEvidence_ref']]
            pepref = match['peptide_ref']
            bseq, seq = db[pepref]
            print(''' <PeptideIdentification score_type="" higher_score_better="false" significance_threshold="0" MZ="%f" RT="%s" >
            <PeptideHit score="%f" sequence="%s" charge="%d" aa_before="%s" aa_after="%s"></PeptideHit>
        </PeptideIdentification>
                ''' % (mz, rt, score, seq, z, enzN, enzC), file=fh)

    print('\n</IdentificationRun>\n</IdXML>',  file=fh)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        mzid2pin(fn)
