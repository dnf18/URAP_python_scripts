from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from pathlib import Path
import json


class Reporter:
    def __init__(self, comparison_json: Path):
        self.comparison_json = Path(comparison_json)

        if not self.comparison_json.exists():
            raise FileNotFoundError(f"Comparison JSON not found: {self.comparison_json}")

        with open(self.comparison_json, "r") as f:
            self.data = json.load(f)

        self.out_pdf = self.comparison_json.parent / "comparison_report.pdf"

    # ------------------------------------------------------------------

    def generate_pdf(self, results=None, histograms=None):
        """
        results: optional additional pass/fail metadata
        histograms: unused but kept for API compatibility
        """

        styles = getSampleStyleSheet()
        story = []

        # Title
        title = Paragraph("<b>MEGAlib Dual-Run Comparison Report</b>", styles['Title'])
        story.append(title)
        story.append(Spacer(1, 18))

        # ------------------------------------------------------------------
        # SUMMARY BLOCK
        # ------------------------------------------------------------------

        ref_mu = self.data["reference"]["mean"]
        ref_sigma = self.data["reference"]["sigma"]
        test_mu = self.data["test"]["mean"]
        test_sigma = self.data["test"]["sigma"]
        sigma_diff = abs(ref_sigma - test_sigma)

        ks_stat = self.data["ks_test"]["statistic"]
        ks_pvalue = self.data["ks_test"]["pvalue"]

        summary_par = Paragraph(
            f"<b>Summary:</b><br/>"
            f"Reference σ: {ref_sigma:.4f}<br/>"
            f"Test σ: {test_sigma:.4f}<br/>"
            f"|Δσ| = {sigma_diff:.4f}<br/><br/>"
            f"KS Test Statistic: {ks_stat:.4f}<br/>"
            f"KS Test p-value: {ks_pvalue:.4f}<br/>",
            styles['BodyText']
        )
        story.append(summary_par)
        story.append(Spacer(1, 12))

        # ------------------------------------------------------------------
        # PASS / FAIL BOX
        # ------------------------------------------------------------------
        
        # A simple rule: pass if p > 0.05 and |Δσ| < 0.1 * reference sigma
        pass_cond = (ks_pvalue > 0.05) and (sigma_diff < 0.10 * ref_sigma)
        pass_text = "PASS" if pass_cond else "FAIL"
        box_color = colors.green if pass_cond else colors.red

        table = Table(
            [[pass_text]],
            colWidths=[80],
            rowHeights=[25]
        )
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), box_color),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.white),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('FONTSIZE', (0, 0), (-1, -1), 14),
            ('BOX', (0, 0), (-1, -1), 1, colors.black),
        ]))
        story.append(table)
        story.append(Spacer(1, 18))

        # ------------------------------------------------------------------
        # HISTOGRAM IMAGE
        # ------------------------------------------------------------------

        img_path = self.data.get("overlay_image", None)
        if img_path and Path(img_path).exists():
            story.append(Paragraph("<b>Energy Histogram Comparison</b>", styles["Heading2"]))
            story.append(Spacer(1, 10))

            img = Image(str(img_path), width=450, height=300)
            story.append(img)
            story.append(Spacer(1, 20))
        else:
            story.append(Paragraph("<b>No overlay image available.</b>", styles["BodyText"]))

        # ------------------------------------------------------------------
        # SAVE PDF
        # ------------------------------------------------------------------

        doc = SimpleDocTemplate(str(self.out_pdf), pagesize=letter)
        doc.build(story)

        print(f"[Reporter] PDF written to: {self.out_pdf}")
        
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

            print("📧 Email sent successfully to:", ", ".join(self.email_recipients))

        except Exception as e:
            print("Failed to send email:", e)

