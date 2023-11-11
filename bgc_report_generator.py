import jinja2

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "bgc_report_template.html.jinja"
template = templateEnv.get_template(TEMPLATE_FILE)


def fill_template(query, genotypes):
    return template.render(query=query, genotypes=genotypes)


def write_report(jinja_output, filename):
    f = open(filename, 'w')
    f.write(jinja_output)
    f.close()
