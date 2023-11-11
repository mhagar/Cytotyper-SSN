from flask import render_template, send_file, flash
from app import app
from app.forms import TestForm
import os.path
import filewrap
import cytotyper_test
import query_parser
import bgc_report_generator


def allowed_upload(filename):
    if '.' in filename:
        if filename.rsplit('.', 1)[1].lower() == 'sqlite':
            return True
    else:
        return False


@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = TestForm()
    rawdata = form.user_upload.data
    rawquery = form.user_upload_query.data
    if form.validate_on_submit() and allowed_upload(rawdata.filename):
        try:
            labels, pfams = query_parser.csv(rawquery.filename)

            qinput = cytotyper_test.query(pfams, labels)

            cytotyper_test.gap_unit = form.user_min_gap_size.data
            cytotyper_test.min_occurence = form.user_min_occurences.data
            analysis = cytotyper_test.analysis(rawdata, qinput)

            report = bgc_report_generator.fill_template(query=analysis.query,
                                                        genotypes=analysis.genotypes)

            bgc_report_generator.write_report(jinja_output=report,
                                              filename='app/output/report.html')

            analysis.cytotable.to_csv('app/output/cytoscape_keys.csv')
            rawdata.save(os.path.join(app.config['UPLOAD_DIR'], rawdata.filename))
            filewrap.rmfile('app/output.zip')
            filewrap.zip_wrap('app/output')
            return send_file('output.zip', attachment_filename='results.zip')
        except Exception as e:
            print(e)
            flash("Either something's wrong with your submission, or there's a bug!"
                  f"Let Mostafa know about this. Tell him it's about '{e}'")
    return render_template('index.html',
                           title='Home',
                           form=form)
