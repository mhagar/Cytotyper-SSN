import cytotyper_test
import bgc_report_generator


class test_input:
    def __init__(self, filename):
        self.filename = filename


finput = test_input('StzF-full.sqlite')

pfams = {1: ['PF06325', 'IPR025799', 'IPR029063'],
         2: ['PF00355', 'IPR036922', 'IPR017941'],
         3: ['PF13535', 'IPR011761'],
         4: ['PF01425', 'IPR036928', 'IPR023631'],
         5: ['PF07859'],
         6: ['PF07690']}


labels = {1: 'PrmA',
          2: 'Reiske',
          3: 'ATP-grasp',
          4: 'Amidase',
          5: 'AB_Hydrolase_3',
          6: 'MFS'}

qinput = cytotyper_test.query(pfams, labels)

cytotyper_test.gap_unit = 500
cytotyper_test.min_occurence = 2
analysis = cytotyper_test.analysis(finput, qinput)

report = bgc_report_generator.fill_template(query=analysis.query,
                                            genotypes=analysis.genotypes)

bgc_report_generator.write_report(jinja_output=report,
                                  filename='test_jan_24.html')
