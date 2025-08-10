from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet

def export_pdf(filename, analysis_data):
    doc = SimpleDocTemplate(filename)
    styles = getSampleStyleSheet()
    elements = []

    for key, value in analysis_data.items():
        elements.append(Paragraph(f"<b>{key}</b>: {value}", styles['Normal']))
        elements.append(Spacer(1, 12))

    doc.build(elements)
    return filename
