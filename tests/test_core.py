import json
from pathlib import Path
import numpy as np
import pysam
import pytest
from scipy.stats import nbinom
from seqgain.io import priority, fragments
from seqgain.kmers import canonical_kmers, canonical_kmer_evidence, analyze_kmers, discovery_forecast
from seqgain.bam import analyze_bam
from seqgain.models import breadth, fit_breadth
from seqgain.cli import main


def test_bam_read_recovery_and_fastq_comparison(tmp_path):
    from seqgain.io import bam_fragments, compare_bam_fastq
    path = tmp_path/'reads.bam'
    with pysam.AlignmentFile(path, 'wb', header={'SQ':[{'SN':'ref','LN':100}]}) as handle:
        for flag, seq in [(81,'GTTT'), (133,'AAAC'), (256,'CCCC')]:
            read=pysam.AlignedSegment()
            read.query_name='pair'
            read.flag=flag
            read.query_sequence=seq
            read.query_qualities=[40]*4
            if not read.is_unmapped:
                read.reference_id=0
                read.reference_start=0
                read.cigarstring='4M'
            handle.write(read)
    a,b=tmp_path/'a.fq',tmp_path/'b.fq'
    fq(a,[('pair/1','AAAC')]);fq(b,[('pair/2','AAAC')])
    recovered=list(bam_fragments(path,tmp_path))
    assert len(recovered)==1
    assert [r[1] for r in recovered[0][1]]==['AAAC','AAAC']
    assert compare_bam_fastq(path,a,b)['match']
    fq(b,[('pair/2','CCCC')])
    comparison=compare_bam_fastq(path,a,b)
    assert comparison['bam']['records']==comparison['fastq']['records']==2
    assert not comparison['match']


def test_bam_missing_qualities_is_explicit(tmp_path):
    from seqgain.io import bam_fragments
    path=tmp_path/'missing.bam'
    with pysam.AlignmentFile(path,'wb',header={'SQ':[{'SN':'ref','LN':100}]}) as handle:
        read=pysam.AlignedSegment();read.query_name='missing';read.flag=4
        read.query_sequence='AAA';handle.write(read)
    with pytest.raises(ValueError,match='lacks name, sequence or qualities'):
        list(bam_fragments(path,tmp_path))


def test_bam_read_group_sampling_keeps_read_names_distinct(tmp_path):
    from seqgain.io import bam_fragments
    path=tmp_path/'read-groups.bam'
    fractions=[.05,.1,.2,.25,.4,.5,.6,.8,1]
    with pysam.AlignmentFile(path,'wb',header={'SQ':[{'SN':'ref','LN':100}],
                                               'RG':[{'ID':'rg1'},{'ID':'rg2'}]}) as handle:
        for i in range(40):
            for rg in ('rg1','rg2'):
                read=pysam.AlignedSegment()
                read.query_name=f'read{i}'
                read.flag=0
                read.reference_id=0
                read.reference_start=i
                read.mapping_quality=60
                read.cigarstring='4M'
                read.query_sequence='ACGT'
                read.query_qualities=[40]*4
                read.set_tag('RG',rg)
                handle.write(read)
    pysam.index(str(path))
    assert priority(('rg1','read/1'),42)==priority(('rg1','read/2'),42)
    assert priority(('rg1','read0'),42)!=priority(('rg1','read1'),42)
    assert priority(('rg1','read0'),42)!=priority(('rg2','read0'),42)
    result=analyze_kmers(None,None,tmp_path,fractions,3,42,3,20,
                         fragment_stream=bam_fragments(path,tmp_path))
    assert result['fragments']==80
    for row in result['curves']:
        expected=sum(priority((rg,f'read{i}'),42,row['replicate'])<row['fraction']
                     for i in range(40) for rg in ('rg1','rg2'))
        assert row['fragments']==expected
    assert 0<next(row['fragments'] for row in result['curves']
                  if row['replicate']==0 and row['fraction']==.1)<80
    coverage=analyze_bam(path,None,fractions,3,42,30,20)
    for row in coverage['targets'][0]['curves']:
        expected=sum(priority((rg,f'read{i}'),42,row['replicate'])<row['fraction']
                     for i in range(40) for rg in ('rg1','rg2'))
        assert row['mean_fragment_depth']==pytest.approx(4*expected/100)


def test_bam_only_default_and_reference_only_escape(tmp_path):
    from seqgain.demo import make_demo
    inputs=tmp_path/'inputs';make_demo(inputs)
    main(['--bam',str(inputs/'sample.bam'),'--outdir',str(tmp_path/'bam-only'),'--replicates','1'])
    result=json.loads((tmp_path/'bam-only/results.json').read_text())
    assert result['kmers']['source']=='bam'
    assert result['kmers']['fragments']==640
    assert result['bam']['summary']['forecast_efforts'][-1]==100
    main(['--bam',str(inputs/'sample.bam'),'--kmer-source','none','--outdir',str(tmp_path/'reference-only'),'--replicates','1'])
    assert json.loads((tmp_path/'reference-only/results.json').read_text())['kmers'] is None


def fq(path, records):
    path.write_text("".join(f"@{name}\n{seq}\n+\n{'I'*len(seq)}\n" for name,seq in records))


def test_pair_sampling_and_quality():
    assert priority("read/1",42) == priority("read/2",42)
    assert priority("read",42) != priority("read",43)
    assert canonical_kmers([("a","AAAC","IIII"),("b","GTTT","IIII")],3,20)=={"AAA","AAC"}
    assert canonical_kmers([("a","AAAC","I!II")],3,20)==set()
    assert canonical_kmers([("a","AAANAAA","IIIIIII")],3,20)=={"AAA"}
    assert canonical_kmer_evidence([("a","AAA","!!!")],3,0)["AAA"] == 1
    q20 = canonical_kmer_evidence([("a","AAA","555")],3,20)
    assert q20["AAA"] == pytest.approx(1-.99**3)
    best = canonical_kmer_evidence([("a","AAANAAA","555NIII")],3,20)
    assert best["AAA"] == pytest.approx(1-.9999**3)


def test_fastq_pair_validation(tmp_path):
    a,b=tmp_path/'a.fq',tmp_path/'b.fq'
    fq(a,[("x/1","AAA")]);fq(b,[("y/2","AAA")])
    with pytest.raises(ValueError,match="matching"):
        list(fragments(a,b))
    a.write_text("@x\nAAA\n+\nII\n")
    with pytest.raises(ValueError,match="Malformed"):
        list(fragments(a))
    a.write_text("@   \nAAA\n+\nIII\n")
    with pytest.raises(ValueError,match="Malformed"):
        list(fragments(a))


def test_kmer_exact_incidence_and_nested_subsampling(tmp_path):
    a=tmp_path/'a.fq';fq(a,[("x","AAAA"),("y","AAAA"),("z","CCCC"),
                            ("u","ACA"),("v","ACA"),("w","ACA")])
    result=analyze_kmers(a,None,tmp_path,[.1,.3,.6,1],3,42,3,20)
    assert result['spectrum']=={1:1,2:1,3:1}
    assert result['distinct']==3 and result['recurrent_distinct']==2
    assert result['phred_error']['expected_false_singletons_upper'] > 0
    for rep in range(3):
        rows=[r for r in result['curves'] if r['replicate']==rep]
        assert [r['distinct'] for r in rows]==sorted(r['distinct'] for r in rows)
        assert rows[-1]['distinct']==3
        supports={"AAA": ["x", "y"], "CCC": ["z"], "ACA": ["u", "v", "w"]}
        for row in rows:
            incidences=[sum(priority(name,42,rep)<row['fraction'] for name in names)
                        for names in supports.values()]
            assert row['singletons']==sum(n==1 for n in incidences)
            assert row['doubletons']==sum(n==2 for n in incidences)
    assert result['phred_error']['expected_false_singletons_upper'] == pytest.approx(1-.9999**3)
    with pytest.raises(ValueError,match="overwrite"):
        analyze_kmers(a,None,tmp_path,[.1,.3,.6,1],3,42,3,20)


def test_duplicate_fastq_names_rejected(tmp_path):
    a=tmp_path/'a.fq';fq(a,[("x","AAAA"),("x","AAAA")])
    with pytest.raises(ValueError,match="Duplicate"):
        analyze_kmers(a,None,tmp_path,[.1,.3,.6,1],1,42,3,20)


def test_incidence_extrapolation():
    m=discovery_forecast(100,1000,100,50,[1,1.5,2])
    y=[p['distinct'] for p in m['points']]
    assert y[0]==1000 and y[0]<y[1]<y[2]<m['asymptote_lower_bound']
    assert y[2]-y[1]<y[1]-y[0]


def make_bam(path):
    header={"HD":{"SO":"coordinate"},"SQ":[{"SN":"seg","LN":30},{"SN":"absent","LN":10}]}
    with pysam.AlignmentFile(path,"wb",header=header) as bam:
        for name,start,cigar,flag,mapq,qual in [
            ('pair',0,'10M',65,60,40),('pair',5,'10M',129,60,40),
            ('deletion',15,'2M3D2M',0,60,40),('skip',15,'2M3N2M',0,60,40),
            ('duplicate',22,'5M',1024,60,40),('secondary',22,'5M',256,60,40),
            ('supplementary',22,'5M',2048,60,40),('lowmap',22,'5M',0,0,40),
            ('unknownmap',22,'5M',0,255,40),('lowbase',22,'5M',0,60,0),
            ('unmapped',-1,'*',4,0,40)]:
            r=pysam.AlignedSegment();r.query_name=name;r.flag=flag
            r.reference_id=0 if start>=0 else -1;r.reference_start=start;r.mapping_quality=mapq
            r.cigarstring=cigar
            length=4 if cigar in ('2M3D2M','2M3N2M') else 10 if cigar=='10M' else 5
            r.query_sequence='A'*length;r.query_qualities=[qual]*length
            bam.write(r)
    pysam.index(str(path))


def test_bam_overlap_cigar_filters_missing_segment(tmp_path):
    path=tmp_path/'sample.bam';make_bam(path)
    r=analyze_bam(path,None,[.1,.3,.6,1],2,42,30,20)
    seg,absent=r['targets']
    full=seg['curves'][-1]
    assert full['breadth']['1']==pytest.approx(19/30)
    assert full['breadth']['2']==pytest.approx(4/30)
    assert full['mean_fragment_depth']==pytest.approx(23/30)
    assert absent['model']['status']=='no_coverage'
    assert absent['model']['forecast_efforts'] == []
    assert seg['thresholds'][:3] == [1,2,3] and seg['thresholds'][-1] == 100
    assert r['summary']['length'] == 40
    summary_full = next(row for row in r['summary']['curves']
                        if row['replicate'] == 0 and row['fraction'] == 1)
    assert summary_full['mean_fragment_depth'] == pytest.approx(23/40)
    assert summary_full['breadth']['1'] == pytest.approx(19/40)
    assert r['coverage_definition']['unit'] == 'fragment'
    forecasts = seg['model']['forecast_breadth']
    assert set(map(int, forecasts)) == set(range(1, 101))
    for depth, values in forecasts.items():
        assert values[0] == pytest.approx(full['breadth'][depth])
        assert all(left <= right + 1e-12 for left,right in zip(values,values[1:]))
    for index in range(len(seg['model']['forecast_efforts'])):
        values = [forecasts[str(depth)][index] for depth in range(1, 101)]
        assert all(left + 1e-12 >= right for left,right in zip(values,values[1:]))
    assert sum(seg['full_depth_histogram'].values())==30
    assert r['record_counts']['primary_records']==9
    supports={'pair':set(range(15)),'deletion':{15,16,20,21},'skip':{15,16,20,21}}
    for row in seg['curves']:
        expected=[sum(pos in positions and priority(name,42,row['replicate']) < row['fraction']
                      for name,positions in supports.items()) for pos in range(30)]
        assert row['mean_fragment_depth']==pytest.approx(sum(expected)/30)
        for depth,value in row['breadth'].items():
            assert value==pytest.approx(sum(d>=int(depth) for d in expected)/30)


def test_pileup_cap_is_not_silent(tmp_path):
    path=tmp_path/'sample.bam';make_bam(path)
    with pytest.raises(ValueError,match="cap reached"):
        analyze_bam(path,None,[.1,.3,.6,1],1,42,30,20,max_depth=1)


def test_target_selection_validation(tmp_path):
    from seqgain.bam import read_targets
    t=tmp_path/'targets.tsv'
    t.write_text('contig\tgroup\tcopies\nmissing\tvirus\t1\n')
    with pytest.raises(ValueError,match='Unknown'):
        read_targets(t,['seg'],[30])
    t.write_text('contig\tgroup\tcopies\nseg\tvirus\t3\n')
    assert read_targets(t,['seg','background'],[30,100])==[
        {'contig':'seg','group':'virus','length':30}]
    t.write_text('contig\nseg\n')
    assert read_targets(t,['seg'],[30]) == [{'contig':'seg','group':'reference','length':30}]


def test_compressed_paired_fastq(tmp_path):
    import gzip
    a,b=tmp_path/'a.fq.gz',tmp_path/'b.fq.gz'
    with gzip.open(a,'wt') as f:
        f.write('@pair/1\nAAAC\n+\nIIII\n')
    with gzip.open(b,'wt') as f:
        f.write('@pair/2\nGTTT\n+\nIIII\n')
    result=analyze_kmers(a,b,tmp_path,[.1,.3,.6,1],1,42,3,20)
    assert result['fragments']==1 and result['bases']==8
    assert result['spectrum']=={1:2}


def test_model_known_scenario_and_holdout():
    rows=[]
    thresholds=[1,5,10,20]
    accessible,mean,shape=.92,22,1.8
    def histogram(fraction):
        probability=shape/(shape+mean*fraction)
        result={depth:accessible*float(nbinom.pmf(depth,shape,probability))
                for depth in range(201)}
        result[0] += 1-accessible
        return result
    for f in [.05,.1,.2,.4,.6,.8,1]:
        rows.append({'replicate':0,'fraction':f,
                     'breadth':{str(d):float(breadth(f,d,accessible,mean,shape)) for d in thresholds}})
    histograms={(0,f):histogram(f) for f in [.05,.1,.2,.4,.6,.8,1]}
    fitted=fit_breadth(rows,thresholds,histograms[(0,1)],histograms,list(range(1,101)))
    assert fitted['accessible_fraction']==pytest.approx(.92,abs=.001)
    assert fitted['holdout_max_absolute_error']<.001
    assert all(0<=v<=1 for values in fitted['forecast_breadth'].values() for v in values)
    assert fitted['forecast_breadth']['10'][fitted['forecast_efforts'].index(2)]==pytest.approx(float(breadth(2,10,.92,22,1.8)),abs=.001)
    assert fitted['retrospective_backtests'][0]['max_absolute_error'] < .001
    assert fitted['forecast_breadth']['10'][0] == pytest.approx(float(breadth(1,10,.92,22,1.8)), abs=.001)


def test_end_to_end_and_html_escaping(tmp_path):
    from seqgain.demo import make_demo
    inputs=tmp_path/'inputs';out=tmp_path/'out';make_demo(inputs)
    main(['--reads1',str(inputs/'reads.fastq'),'--bam',str(inputs/'sample.bam'),
          '--outdir',str(out),'--sample','</script><script>alert(1)</script>','--replicates','1'])
    result=json.loads((out/'results.json').read_text())
    assert result['kmers']['fragments']==640
    assert result['kmers']['source']=='bam'
    assert result['provenance']['bam_fastq_comparison']['match']
    assert result['kmers']['forecast']['points'][-1]['effort']==100
    assert result['schema_version']==4
    assert 'required_effort' not in result['bam']['targets'][0]['model']
    assert 'target_depth' not in result['bam']['targets'][0]
    assert result['kmers']['retrospective_backtests']
    assert len(result['bam']['targets'])==2
    assert result['bam']['summary']['label']=='All references'
    report=(out/'report.html').read_text()
    assert '</script><script>alert(1)</script>' not in report
    # The bundled library contains dormant CDN defaults; no external scripts are loaded.
    from html.parser import HTMLParser
    class Scripts(HTMLParser):
        sources = []
        def handle_starttag(self, tag, attrs):
            if tag == 'script':
                self.sources.extend(v for k,v in attrs if k == 'src')
    parsed = Scripts(); parsed.feed(report)
    assert parsed.sources == []
    assert 'type="range" min="1" max="100"' in report
    assert 'All references' in report
    assert 'Minimum depth' in report
    assert 'REFERENCE TARGET SCENARIO' not in report
    assert "'Anchored '+n+'× scenario'" in report
    assert not (out/'kmers.sqlite').exists()
    assert (out/'curves.tsv').stat().st_size>100
    check_report_script(out/'report.html')


@pytest.mark.parametrize('extra', [['--fractions','.1,.2,.3'],['--replicates','0'],['--baseq','100'],['--explore-filters'],['--kmer-source','none']])
def test_cli_validation(tmp_path,extra):
    a=tmp_path/'a.fq';fq(a,[('x','AAAA')])
    with pytest.raises(SystemExit) as exc:
        main(['analyze','--reads1',str(a),'--outdir',str(tmp_path/'out'),*extra])
    assert exc.value.code==2


def test_progress_heartbeat_and_cleanup(monkeypatch):
    import io
    import threading
    from seqgain.progress import Progress
    monkeypatch.setenv('COLUMNS', '240')

    class Output(io.StringIO):
        def isatty(self):
            return True

        def write(self, value):
            result = super().write(value)
            if "ETA ~" in value:
                emitted.set()
            return result

    emitted = threading.Event()
    output = Output()
    reporter = Progress(enabled=True, interval=.01, stream=output)
    with reporter:
        reporter.phase('Reads', 10, 'reads')
        reporter.advance(5)
        assert emitted.wait(2), 'heartbeat must report while computation is busy'
    assert '5/10 reads (50.0%)' in output.getvalue()
    assert 'elapsed' in output.getvalue() and 'reads/s' in output.getvalue()
    assert not reporter._thread.is_alive()
    assert output.getvalue().count('\n') == 1
    assert '\r\033[2K' in output.getvalue()
    output = Output()
    with pytest.raises(ValueError):
        with Progress(enabled=True, stream=output) as reporter:
            reporter.phase('Unknown input', unit='reads')
            reporter.advance(3)
            raise ValueError('failure')
    assert 'ETA unavailable' in output.getvalue()
    assert 'stopped' in output.getvalue()
    assert '(100.0%)' not in output.getvalue()
    assert not reporter._thread.is_alive()


def test_bam_progress_and_quiet(tmp_path, capsys):
    path=tmp_path/'sample.bam';make_bam(path)
    main(['--bam',str(path),'--outdir',str(tmp_path/'progress'),'--k','3','--replicates','1'])
    stderr=capsys.readouterr().err
    assert '11/11 records (100.0%)' in stderr
    assert '9/9 reads (100.0%)' in stderr
    assert '30/30 positions (100.0%)' in stderr
    assert 'Fitting coverage forecast' in stderr
    main(['--bam',str(path),'--outdir',str(tmp_path/'quiet'),'--kmer-source','none',
          '--replicates','1','--quiet'])
    assert capsys.readouterr().err == ''


def test_fastq_only_progress_and_report(tmp_path, capsys):
    reads=tmp_path/'reads.fq';fq(reads,[('one','AAAC'),('two','AAAA')])
    main(['--reads1',str(reads),'--outdir',str(tmp_path/'fastq'),
          '--k','3','--replicates','1'])
    progress=capsys.readouterr().err
    assert '2 reads' in progress and 'done' in progress
    check_report_script(tmp_path/'fastq/report.html')


def test_progress_log_has_one_line_per_phase():
    import io
    from seqgain.progress import Progress
    output = io.StringIO()
    with Progress(enabled=True, stream=output) as reporter:
        reporter.phase('Reads', 10, 'reads')
        reporter.advance(5)
        reporter._emit()
        reporter.advance(5)
        reporter.phase('Fitting')
        reporter._emit()
        reporter.finish()
        reporter.finish()
    lines = output.getvalue().splitlines()
    assert len(lines) == 2
    assert '10/10 reads (100.0%)' in lines[0]
    assert '[Fitting]' in lines[1]
    assert '\r' not in output.getvalue() and '\033' not in output.getvalue()
    assert reporter._thread is None


def test_progress_terminal_width_and_phase_lines(monkeypatch):
    import io
    from seqgain.progress import Progress
    class Terminal(io.StringIO):
        def isatty(self):
            return True
    monkeypatch.setenv('COLUMNS', '60')
    output = Terminal()
    with Progress(enabled=True, stream=output) as reporter:
        reporter.phase('Coverage ' + 'long_contig_name'*20, 100, 'positions')
        reporter.update(50)
        reporter._emit()
        reporter.phase('Fitting')
        reporter.finish()
        output.write('report.html\n')
    assert output.getvalue().count('\n') == 3
    assert '50.0%' in output.getvalue()
    for update in output.getvalue().split('\r\033[2K')[1:]:
        assert len(update.split('\n')[0]) <= 59


def test_explorable_filters(tmp_path):
    path=tmp_path/'sample.bam';make_bam(path)
    main(['--bam',str(path),'--outdir',str(tmp_path/'filters'),'--kmer-source','none',
          '--replicates','1','--fractions','.1,.3,.6,1','--mapq','25','--include-duplicates',
          '--explore-filters','--quiet'])
    result=json.loads((tmp_path/'filters/results.json').read_text())
    variants=[result['bam'],*result['bam_filter_variants']]
    assert len(variants)==8
    assert result['bam']['coverage_definition']['minimum_mapq']==25
    assert result['bam']['coverage_definition']['include_marked_duplicates']
    assert {(v['coverage_definition']['minimum_mapq'],v['coverage_definition']['include_marked_duplicates'])
            for v in variants} == {(q,dup) for q in [0,20,25,30] for dup in [False,True]}
    for variant in variants:
        definition=variant['coverage_definition']
        mapq=definition['minimum_mapq'];duplicates=definition['include_marked_duplicates']
        expected=23+5*(mapq==0)+5*duplicates
        full=variant['summary']['curves'][-1]
        assert full['mean_fragment_depth']==pytest.approx(expected/40)
        assert variant['record_counts']['filtered_unknown_mapq_records']==1
        assert variant['targets'][1]['model']['status']=='no_coverage'
        segment=variant['targets'][0]
        for depth,values in segment['model']['forecast_breadth'].items():
            assert values[0]==pytest.approx(segment['curves'][-1]['breadth'][depth])
    default=next(v for v in variants if v['coverage_definition']['minimum_mapq']==30
                 and not v['coverage_definition']['include_marked_duplicates'])
    standalone=analyze_bam(path,None,[.1,.3,.6,1],1,42,30,20)
    assert default==json.loads(json.dumps(standalone))
    check_report_script(tmp_path/'filters/report.html')


def check_report_script(path):
    """Optional frontend regression check; Node is not an analysis dependency."""
    import shutil
    import subprocess
    node = shutil.which('node')
    if node:
        subprocess.run([node, str(Path(__file__).with_name('check_report.cjs')), str(path)],
                       check=True, capture_output=True, text=True, timeout=30)
