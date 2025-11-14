#!/usr/bin/env python3

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


class Reporter:

    #  Inputs validation results and histograms from Comparator (not written yet)
    #  Outputs PDF summary report, send report through email (optional)
    #  Supervisor takes results from comparator and calls reporter for outputs

    def __init__(self, config_json: str = None, output_dir: str = "./reports", email_recipients=None):

        # Initialize Reporter 
    
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load JSON if provided in Steering
        self.config = {}
        if config_json and os.path.exists(config_json):
            with open(config_json, "r") as f:
                self.config = json.load(f)

        self.email_recipients = email_recipients or []

    def generate_pdf(self, results: dict, histograms: dict = None, output_name: str = "Validation_Report.pdf"):
        
        # Creates PDF report (assumes results and histograms dict from Comparator for now)
        
        output_path = self.output_dir / output_name
        doc = SimpleDocTemplate(str(output_path), pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        # Title
        story.append(Paragraph("<b>MEGAlib Validation Report</b>", styles["Title"]))
        story.append(Spacer(1, 12))

        # Configuration Info 
        story.append(Paragraph("<b>Configuration Summary</b>", styles["Heading2"]))
        config_data = [[k, str(v)] for k, v in self.config.items()]
        config_table = Table(config_data, colWidths=[150, 350])
        config_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ]))
        story.append(config_table)
        story.append(Spacer(1, 12))

        # Results Summary 
        story.append(Paragraph("<b>Validation Results</b>", styles["Heading2"]))
        result_data = [[k, str(v)] for k, v in results.items()]
        result_table = Table(result_data, colWidths=[200, 300])
        result_table.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ]))
        story.append(result_table)
        story.append(Spacer(1, 12))

        # Histograms (embedded plots) 
        if histograms:
            story.append(Paragraph("<b>Comparison Histograms</b>", styles["Heading2"]))
            for label, img_path in histograms.items():
                if os.path.exists(img_path):
                    story.append(Paragraph(label, styles["Normal"]))
                    story.append(Image(img_path, width=400, height=250))
                    story.append(Spacer(1, 12))

        # Final Status
        status = "PASS" if results.get("pass") else "FAIL"
        story.append(Paragraph(f"<b>Overall Test Status: {status}</b>", styles["Heading2"]))

        # Build PDF
        doc.build(story)
        print(f"PDF report generated at: {output_path}")
        return output_path

    def send_email(self, subject: str, body: str, attachment_path: str, sender_email: str, sender_password: str):
        
        # Emails the report to recipients
        
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
            msg.add_attachment(f.read(), maintype="application", subtype="pdf", filename=os.path.basename(attachment_path))

        # Send via SMTP
        try:
            with smtplib.SMTP_SSL("smtp.gmail.com", 465) as smtp:
                smtp.login(sender_email, sender_password)
                smtp.send_message(msg)
            print("📧 Email sent successfully to:", ", ".join(self.email_recipients))
        except Exception as e:
            print("Failed to send email:", e)