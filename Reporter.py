import os
import json
import smtplib
import matplotlib.pyplot as plt
from email.message import EmailMessage
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors

"""""
This is the Reporter Class

This is where we attatch all relevant data/graphs onto a PDF aptly named Validation_Report.pdf
"""""

class Reporter:

    # Inputs validation results and histograms from Comparator
    # Outputs PDF summary report, send report through email (optional)

    def __init__(self, config_json: str = None, output_dir: str = "./reports", email_recipients=None):

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load JSON if provided in Steering
        self.config = {}
        if config_json and os.path.exists(config_json):
            with open(config_json, "r") as f:
                self.config = json.load(f)

        self.email_recipients = email_recipients or []

    def generate_pdf(self, results: dict, histograms: dict = None, output_name: str = "Validation_Report.pdf"):
        """
        Creates PDF report (uses results and histogram paths from Comparator).
        """

        output_path = self.output_dir / output_name
        doc = SimpleDocTemplate(str(output_path), pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        # Title

        story.append(Paragraph("<b>MEGAlib Validation Report</b>", styles["Title"]))
        story.append(Spacer(1, 14))

        # Configuration Table

        story.append(Paragraph("<b>Configuration Summary</b>", styles["Heading2"]))

        if len(self.config) > 0:
            config_data = [["Parameter", "Value"]] + [[k, str(v)] for k, v in self.config.items()]
            config_table = Table(config_data, colWidths=[150, 350])
            config_table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
            ]))
            story.append(config_table)
            story.append(Spacer(1, 12))
        else:
            story.append(Paragraph("No configuration provided.", styles["BodyText"]))
            story.append(Spacer(1, 12))

        # Results Table

        story.append(Paragraph("<b>Validation Results</b>", styles["Heading2"]))

        result_data = [["Metric", "Value"]]
        for k, v in results.items():
            result_data.append([str(k), str(v)])

        result_table = Table(result_data, colWidths=[200, 300])
        result_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ]))
        story.append(result_table)
        story.append(Spacer(1, 16))

        # Histogram Embedding (side-by-side)

        if histograms:
            story.append(Paragraph("<b>Comparison Histograms</b>", styles["Heading2"]))

            rows = []
            labels = list(histograms.keys())
            paths  = list(histograms.values())

            images = []
            for p in paths:
                if os.path.exists(p):
                    images.append(Image(p, width=240, height=180))
                else:
                    images.append(Paragraph("Missing image", styles["BodyText"]))

            # Pair images two per row
            for i in range(0, len(images), 2):
                row = images[i:i+2]
                if len(row) < 2:
                    row.append("")
                rows.append(row)

            table = Table(rows, colWidths=[260, 260])
            table.setStyle(TableStyle([
                ("GRID", (0, 0), (-1, -1), 0.25, colors.black),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ]))

            story.append(table)
            story.append(Spacer(1, 16))

        # Final PASS/FAIL Status
    
        status = "PASS" if results.get("pass") else "FAIL"
        color = "green" if status == "PASS" else "red"
        story.append(Paragraph(f"<b>Overall Test Status:</b> <font color='{color}'>{status}</font>", styles["Heading2"]))

        # Build the PDF
        doc.build(story)
        print(f"PDF report generated at: {output_path}")
        return output_path

    # Email Sending (original structure preserved)

    def send_email(self, subject: str, body: str, attachment_path: str,
                   sender_email: str, sender_password: str):
        """
        Emails the PDF report to all recipients.
        """

        if not self.email_recipients:
            print("No email recipients specified. Skipping email step.")
            return

        msg = EmailMessage()
        msg["Subject"] = subject
        msg["From"] = sender_email
        msg["To"] = ", ".join(self.email_recipients)
        msg.set_content(body)

        # Attach the PDF
        with open(attachment_path, "rb") as f:
            msg.add_attachment(
                f.read(),
                maintype="application",
                subtype="pdf",
                filename=os.path.basename(attachment_path)
            )

        # Send via SMTP
        try:
            with smtplib.SMTP_SSL("smtp.gmail.com", 465) as smtp:
                smtp.login(sender_email, sender_password)
                smtp.send_message(msg)

            print("Email sent successfully to:", ", ".join(self.email_recipients))

        except Exception as e:
            print("Failed to send email:", e)
