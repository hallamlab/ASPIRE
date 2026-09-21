"""Execute the production cohort-routing block against small Nextflow fixtures."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')


@unittest.skipUnless(NEXTFLOW, 'Nextflow required')
class CohortDependencyTests(unittest.TestCase):
    def run_case(self, produce=True, enabled=True):
        source = (ROOT/'asv_pipeline.nf').read_text()
        start = source.index('    if( diversityEnabled ) {\n        def cohortInput')
        end = source.index('    if( clustermapsEnabled ) {', start)
        routing = source[start:end]
        # Reuse the actual emitted-output declaration from the producer.
        module = (ROOT/'workflow/modules/indicators_diagnostics.nf').read_text()
        output = next(line for line in module.splitlines() if 'emit: sample_cohort' in line)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            producer = 'printf "sampleID\\ns1\\n" > sample_analysis_cohort.tsv' if produce else 'true'
            nf = '''nextflow.enable.dsl=2
process COMMUNITY_PREDICTOR_COMPARISON {
 output:
OUTPUT
 script:
 """
PRODUCER
 """
}
process DIVERSITY_ANALYSIS {
 input:
 val(metadata)
 val(counts)
 path(cohort)
 output:
 stdout
 script:
 """
 test -s "${cohort}"
 echo MATCHED_COHORT_READY
 """
}
workflow {
 def diversityEnabled = true
 def communityPredictorEnabled = ENABLED
 def diversityMatchedCohortSource = 'community_predictor_comparison'
 def diversityMatchedCohortPath = ''
 def metaMicroForDiversity = 'metadata'
 def asvFinalForDiversity = 'counts'
 if(communityPredictorEnabled) { COMMUNITY_PREDICTOR_COMPARISON() }
ROUTING
 DIVERSITY_ANALYSIS.out.view()
}
'''.replace('OUTPUT',output).replace('PRODUCER',producer).replace('ENABLED',str(enabled).lower()).replace('ROUTING',routing)
            (p/'main.nf').write_text(nf)
            result = subprocess.run([NEXTFLOW,'run','main.nf','-ansi-log','false'],cwd=p,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
            return result

    def test_fresh_generated_cohort(self):
        result=self.run_case()
        self.assertEqual(result.returncode,0,result.stdout)
        self.assertIn('MATCHED_COHORT_READY',result.stdout)

    def test_missing_cohort_fails_explicitly(self):
        result=self.run_case(produce=False)
        self.assertNotEqual(result.returncode,0,result.stdout)
        self.assertIn('produced no sample_analysis_cohort.tsv',result.stdout)

    def test_disabled_producer_fails_explicitly(self):
        result=self.run_case(enabled=False)
        self.assertNotEqual(result.returncode,0,result.stdout)
        self.assertIn('requires community_predictor_comparison.enabled=true',result.stdout)
